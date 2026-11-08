
        
# --- Binary file helpers ---

"""
    _is_binary(path)

Return `true` if `path` has a binary extension (`.bin` or `.bin.gz`).
"""
_is_binary(path::AbstractString)::Bool = endswith(path, ".bin") || endswith(path, ".bin.gz")

"""
    _open_output(path)

Open an output file for reading, transparently handling gzip compression.
"""
_open_output(path::AbstractString)::IO = endswith(path, ".gz") ? GZip.open(path) : open(path)

"""
    _read_bin_header(io)

Read and validate the 16-byte binary header from an IO stream.
Returns the number of Float64 values per row (`ncols`).
"""
function _read_bin_header(io::IO)::Int
    magic = read(io, UInt32)
    magic == GRN._BIN_MAGIC || error(
        "Not a GRN binary file (bad magic: 0x$(string(magic, base=16, pad=8)))")
    version = read(io, UInt32)
    version == GRN._BIN_VERSION || error(
        "Unsupported GRN binary version: $version (expected $(GRN._BIN_VERSION))")
    ncols = read(io, Int64)
    ncols > 0 || error("Invalid GRN binary file: ncols = $ncols")
    return ncols
end

"""
    read_binary(path)

Read a GRN binary output file (written with `format=:bin`) and return a
matrix with one row per MCMC iteration and one column per parameter value.
Handles both plain `.bin` and gzip-compressed `.bin.gz` files.

# Example
```julia
# After running grn with format=:bin:
G = GRN.read_binary("output/G_effects.bin")       # 51 × 25 matrix
Q = GRN.read_binary("output/Q_effects.bin.gz")    # 51 × 50 matrix
```
"""
function read_binary(path::AbstractString)
    io = _open_output(path)
    try
        ncols = _read_bin_header(io)
        data_bytes = read(io)
        n_bytes = length(data_bytes)
        n_bytes % sizeof(Float64) == 0 || error(
            "Corrupt binary file: data size ($n_bytes bytes) is not a multiple of 8")
        nvals = n_bytes ÷ sizeof(Float64)
        nvals % ncols == 0 || error(
            "Corrupt binary file: $nvals values is not a multiple of ncols=$ncols")
        nrows = nvals ÷ ncols
        return reshape(reinterpret(Float64, data_bytes), ncols, nrows)'
    finally
        close(io)
    end
end

"""
    _posterior_mean_bin(path)

Streaming posterior mean for binary files. Reads row by row without
loading the entire file into memory.
"""
function _posterior_mean_bin(path::AbstractString)::Vector{Float64}
    io = _open_output(path)
    try
        ncols = _read_bin_header(io)
        sum = zeros(Float64, ncols)
        row = Vector{Float64}(undef, ncols)
        nrows = 0
        while !eof(io)
            read!(io, row)
            @inbounds for i in 1:ncols
                sum[i] += row[i]
            end
            nrows += 1
        end
        nrows == 0 && error("No data rows found in $path")
        return sum ./ nrows
    finally
        close(io)
    end
end

"""
    _pip_bin(path)

Streaming posterior inclusion probability for binary files.
"""
function _pip_bin(path::AbstractString)::Vector{Float64}
    io = _open_output(path)
    try
        ncols = _read_bin_header(io)
        counts = zeros(Int, ncols)
        row = Vector{Float64}(undef, ncols)
        nrows = 0
        while !eof(io)
            read!(io, row)
            @inbounds for i in 1:ncols
                row[i] != 0.0 && (counts[i] += 1)
            end
            nrows += 1
        end
        nrows == 0 && error("No data rows found in $path")
        return counts ./ nrows
    finally
        close(io)
    end
end

# --- Readers for files that are still being written ---
#
# While `grn` runs, output files end wherever the last 1 MB buffer flush
# stopped: the last row may lack its end-of-line character (and may be cut
# mid-number), a binary row may be incomplete, a binary header may not be on
# disk yet, and a gzip stream has no trailer. GZip.jl throws on a gzip stream
# without trailer and discards the bytes it has read, so gzip files are read
# with zlib's own `gzread`, which returns the available data and then 0.

"""
    _PartialGz(file)

Handle to a gzip file opened with zlib's `gzopen` for reading files that may
still be in progress. Read with `_read_partial!`, close with `_close_partial`.
"""
struct _PartialGz
    file::Ptr{Cvoid}
end

"""
    _open_partial(path)

Open an output file for reading while it may still be written to. Gzip files
(`.gz`) are opened with zlib directly; other files with `open`.
"""
function _open_partial(path::AbstractString)::Union{IOStream, _PartialGz}
    endswith(path, ".gz") || return open(path, "r")
    isfile(path) || throw(SystemError("opening file $path", 2))
    f = @static if Sys.iswindows()
        ccall((:gzopen_w, libz), Ptr{Cvoid}, (Cwstring, Cstring), path, "rb")
    else
        ccall((:gzopen, libz), Ptr{Cvoid}, (Cstring, Cstring), path, "rb")
    end
    f == C_NULL && error("Could not open gzip file $path")
    return _PartialGz(f)
end

"""
    _read_partial!(src, buf, nb=length(buf))

Read `nb` bytes into the start of `buf`, or fewer if the end of the available
data is reached. Returns the number of bytes read; fewer than `nb` means the
end of the data written so far. A gzip stream that ends without a trailer is
treated as the end of the data, not as an error.
"""
function _read_partial!(src::IOStream, buf::Vector{UInt8}, nb::Int=length(buf))::Int
    return readbytes!(src, buf, nb)
end

function _read_partial!(src::_PartialGz, buf::Vector{UInt8}, nb::Int=length(buf))::Int
    nb <= length(buf) || throw(ArgumentError("nb = $nb exceeds the buffer length"))
    nr = 0
    while nr < nb
        # gzread returns a C int, so read at most 1 GiB per call
        len = min(nb - nr, 1 << 30)
        ret = ccall((:gzread, libz), Cint, (Ptr{Cvoid}, Ptr{UInt8}, Cuint),
            src.file, pointer(buf) + nr, len)
        ret < 0 && error("Corrupt gzip data (zlib gzread returned $ret)")
        ret == 0 && break
        nr += ret
    end
    return nr
end

"""
    _close_partial(src)

Close a source opened with `_open_partial`. The zlib status for a gzip stream
without trailer is ignored, since an unfinished file is expected here.
"""
_close_partial(src::IOStream)::Nothing = (close(src); nothing)
_close_partial(src::_PartialGz)::Nothing =
    (ccall((:gzclose, libz), Cint, (Ptr{Cvoid},), src.file); nothing)

"""
    _partial_clean(src)

After reading to the end of the available data, return `true` if the source
ended cleanly. A plain file always does. For a gzip file, zlib reports
`Z_OK` when the last gzip stream was finished (trailer present) and
`Z_BUF_ERROR` when the data ended inside a stream, i.e. the file is still
being written or the writing process was stopped before closing it.
"""
_partial_clean(::IOStream)::Bool = true
function _partial_clean(src::_PartialGz)::Bool
    err = Ref{Cint}(0)
    ccall((:gzerror, libz), Cstring, (Ptr{Cvoid}, Ref{Cint}), src.file, err)
    return err[] == 0
end

"""
    _accumulate!(sums, counts, row)

Add one complete row to the running sums and non-zero counts.
"""
function _accumulate!(sums::Vector{Float64}, counts::Vector{Int},
    row::AbstractVector{Float64})::Nothing
    @inbounds for i in eachindex(sums)
        v = row[i]
        sums[i] += v
        v != 0.0 && (counts[i] += 1)
    end
    return nothing
end

"""
    _preliminary_result(sums, counts, nrows)

Turn running sums and non-zero counts into the NamedTuple returned by
`preliminary_posterior`.
"""
function _preliminary_result(sums::Vector{Float64}, counts::Vector{Int},
    nrows::Int)::NamedTuple
    nrows == 0 && return (mean=Float64[], pip=Float64[], iterations=0)
    return (mean=sums ./ nrows, pip=counts ./ nrows, iterations=nrows)
end

"""
    _skip_partial!(src, nb, scratch)

Read and discard `nb` bytes through `scratch`. Returns `false` if the
available data ended first.
"""
function _skip_partial!(src::Union{IOStream, _PartialGz}, nb::Int,
    scratch::Vector{UInt8})::Bool
    while nb > 0
        m = min(nb, length(scratch))
        _read_partial!(src, scratch, m) == m || return false
        nb -= m
    end
    return true
end

"""
    _check_columns(cols, ncols, path)

Throw an `ArgumentError` if a requested column (sorted, unique) is beyond the
`ncols` values per row of the file.
"""
function _check_columns(cols::Union{Nothing, Vector{Int}}, ncols::Int,
    path::AbstractString)::Nothing
    (cols === nothing || isempty(cols) || cols[end] <= ncols) || throw(ArgumentError(
        "Column $(cols[end]) requested, but $path has $ncols values per row"))
    return nothing
end

"""
    _scan_bin_rows(onrow!, path, cols; complete=nothing)

Call `onrow!(row)` for each complete row of a binary output file that may
still be written, and return the number of rows. With `cols === nothing`,
`row` holds all values of the row; otherwise it holds the values of the
sorted, unique columns `cols`, and the other values are skipped without
being stored. A missing or incomplete header gives zero rows, and an
incomplete last row is ignored. `row` is reused between calls.

With `complete` given (a `Ref{Bool}`, only for `cols === nothing`), it is
set to `true` if the file ended cleanly: no incomplete last row and, for
gzip, a finished stream (see `_partial_clean`).
"""
function _scan_bin_rows(onrow!::F, path::AbstractString,
    cols::Union{Nothing, Vector{Int}};
    complete::Union{Nothing, Base.RefValue{Bool}}=nothing)::Int where {F}
    (complete === nothing || cols === nothing) || throw(ArgumentError(
        "complete is only determined when whole rows are read (cols === nothing)"))
    complete === nothing || (complete[] = false)
    src = _open_partial(path)
    try
        hdr = Vector{UInt8}(undef, _BIN_HEADER_SIZE)
        _read_partial!(src, hdr) < _BIN_HEADER_SIZE && return 0
        ncols = _read_bin_header(IOBuffer(hdr))
        _check_columns(cols, ncols, path)
        nrows = 0

        if cols === nothing
            # Whole rows
            rowbytes = Vector{UInt8}(undef, ncols * sizeof(Float64))
            row = reinterpret(Float64, rowbytes)
            nlast = 0  # bytes read for the row after the last complete one
            while (nlast = _read_partial!(src, rowbytes)) == length(rowbytes)
                onrow!(row)
                nrows += 1
            end
            complete === nothing || (complete[] = nlast == 0 && _partial_clean(src))
        else
            # Selected columns: read 8 bytes per wanted value, skip the rest
            row = zeros(Float64, length(cols))
            val = zeros(Float64, 1)
            valbytes = Vector{UInt8}(undef, sizeof(Float64))
            scratch = Vector{UInt8}(undef, min(_BUFSIZE, ncols * sizeof(Float64)))
            while true
                pos = 0  # values of the current row consumed so far
                complete = true
                for k in eachindex(cols)
                    j = cols[k]
                    if !_skip_partial!(src, (j - 1 - pos) * sizeof(Float64), scratch) ||
                       _read_partial!(src, valbytes) < sizeof(Float64)
                        complete = false
                        break
                    end
                    copyto!(val, reinterpret(Float64, valbytes))
                    row[k] = val[1]
                    pos = j
                end
                complete = complete &&
                    _skip_partial!(src, (ncols - pos) * sizeof(Float64), scratch)
                complete || break
                onrow!(row)
                nrows += 1
            end
        end
        return nrows
    finally
        _close_partial(src)
    end
end

"""
    _scan_tsv_rows(onrow!, path, delim, cols; complete=nothing)

Call `onrow!(row)` for each complete row of a text output file that may
still be written, and return the number of rows. Only rows terminated by an
end-of-line character are used; an unterminated last row (possibly cut
mid-number) is ignored. With `cols === nothing`, `row` holds all values of
the row; otherwise it holds the values of the sorted, unique columns `cols`,
and only those fields are parsed. The file is read in 1 MB chunks, so memory
use is one row of values, not one row of text. Blank lines and `\\r\\n` line
endings are accepted. `row` is reused between calls.

With `complete` given (a `Ref{Bool}`), it is set to `true` if the file ended
cleanly: no unterminated last row and, for gzip, a finished stream (see
`_partial_clean`).
"""
function _scan_tsv_rows(onrow!::F, path::AbstractString, delim::Char,
    cols::Union{Nothing, Vector{Int}};
    complete::Union{Nothing, Base.RefValue{Bool}}=nothing)::Int where {F}
    isascii(delim) && delim != '\n' || throw(ArgumentError(
        "delim must be an ASCII character other than newline, got $(repr(delim))"))
    d = UInt8(delim)
    nl = UInt8('\n')

    src = _open_partial(path)
    chunk = Vector{UInt8}(undef, _BUFSIZE)
    carry = UInt8[]     # bytes of a field that continues in the next chunk
    row = cols === nothing ? Float64[] : zeros(Float64, length(cols))
    ncols = 0           # values per row, set by the first complete row
    col = 0             # values read so far in the current row
    k = 1               # next wanted entry of `cols`
    nrows = 0

    try
        while true
            n = _read_partial!(src, chunk)
            n == 0 && break

            # One String per chunk, so fields are parsed as non-copying SubStrings
            s = String(chunk[1:n])
            start = 1
            for i in 1:n
                b = @inbounds chunk[i]
                (b == d || b == nl) || continue

                # The field is the carried bytes plus bytes start:i-1 of this chunk
                if isempty(carry)
                    field = SubString(s, start, i - 1)
                else
                    append!(carry, codeunits(SubString(s, start, i - 1)))
                    field = String(carry)  # empties carry
                end
                start = i + 1

                # Skip blank lines
                b == nl && col == 0 && all(isspace, field) && continue

                col += 1
                nrows > 0 && col > ncols && error(
                    "Row $(nrows + 1) of $path has more than $ncols values")

                # Parse the field if it is wanted
                if cols === nothing || (k <= length(cols) && @inbounds(cols[k]) == col)
                    v = tryparse(Float64, field)
                    v === nothing && error(
                        "Non-numeric value $(repr(field)) in row $(nrows + 1) of $path")
                    if cols === nothing
                        nrows == 0 ? push!(row, v) : (@inbounds row[col] = v)
                    else
                        @inbounds row[k] = v
                        k += 1
                    end
                end

                # End of line: the row is complete
                if b == nl
                    if nrows == 0
                        ncols = col
                        _check_columns(cols, ncols, path)
                    else
                        col == ncols || error(
                            "Row $(nrows + 1) of $path has $col values, expected $ncols")
                    end
                    onrow!(row)
                    nrows += 1
                    col = 0
                    k = 1
                end
            end
            append!(carry, view(chunk, start:n))
        end
        complete === nothing ||
            (complete[] = isempty(carry) && col == 0 && _partial_clean(src))
    finally
        _close_partial(src)
    end

    # Values in `row` and bytes in `carry` belong to an unterminated last row
    # and are dropped
    return nrows
end

"""
    _scan_rows(onrow!, path, delim, cols; complete=nothing)

Dispatch to `_scan_bin_rows` or `_scan_tsv_rows` by file extension.
"""
function _scan_rows(onrow!::F, path::AbstractString, delim::Char,
    cols::Union{Nothing, Vector{Int}};
    complete::Union{Nothing, Base.RefValue{Bool}}=nothing)::Int where {F}
    return _is_binary(path) ? _scan_bin_rows(onrow!, path, cols; complete=complete) :
        _scan_tsv_rows(onrow!, path, delim, cols; complete=complete)
end

"""
    preliminary_posterior(path; delim='\\t')

Posterior mean and posterior inclusion probability from an output file of an
analysis that may still be running. Use it to look at results before `grn`
has finished, or at the files of a run that was stopped.

Only complete iterations are used:
 - text files (`.tsv`, `.tsv.gz`): rows without an end-of-line character
   (typically a row cut off mid-write) are ignored
 - binary files (`.bin`, `.bin.gz`): an incomplete last row is ignored, and a
   header that is not yet on disk gives zero iterations
 - gzip files without a trailer (still being written) are read up to the
   last data that zlib can decompress

The files on disk lag behind the sampler: output passes through a 1 MB
buffer per file and then the stream's own buffer (Julia's `IOStream` buffer,
or zlib's buffer and compressor state for `.gz` files). Small output files,
such as `base_residual_variance` or `G_pi`, can therefore show few or no
iterations long after the sampler has produced them.

The mean and PIP are computed in one pass over the file. For finished files
the result equals `posterior_mean` and `posterior_inclusion_probability`.
Memory use is about one row of values plus two vectors of the same length.

# Returns
A NamedTuple `(mean, pip, iterations)`:
 - `mean`:       Posterior mean of each value (`Vector{Float64}`)
 - `pip`:        Proportion of iterations where each value is non-zero
                 (`Vector{Float64}`); the PIP for effects with a spike-and-slab prior
 - `iterations`: Number of complete iterations used (`Int`). When it is 0,
                 `mean` and `pip` are empty.

# Example
```julia
# While grn(...; outpath="output") is running in another Julia session:
res = preliminary_posterior("output/G_effects.tsv.gz")
res.iterations                       # number of iterations written so far
G_mean = reshape(res.mean, m, m)
G_pip = reshape(res.pip, m, m)
```
"""
function preliminary_posterior(path::AbstractString; delim::Char='\t')::NamedTuple
    # Accumulators are sized on the first row; they are filled in place (not
    # reassigned) so the closure does not box them
    sums = Float64[]
    counts = Int[]
    nrows = _scan_rows(path, delim, nothing) do row
        if isempty(sums)
            append!(sums, zeros(Float64, length(row)))
            append!(counts, zeros(Int, length(row)))
        end
        _accumulate!(sums, counts, row)
    end
    return _preliminary_result(sums, counts, nrows)
end

"""
    read_output(path; delim='\\t')

Read all complete iterations of an output file into a matrix with one row
per iteration and one column per value. Works for `.tsv`, `.tsv.gz`, `.bin`,
and `.bin.gz`, for finished files and for files of an analysis that is
still running or was stopped (incomplete last rows are ignored, as in
`preliminary_posterior`). Returns a 0 × 0 matrix when no iteration is
complete yet.

Meant for the small output files. The whole file is held in memory, so do
not use it on effects files at large scale (one row of `G_effects` is
`m * m` values); use `preliminary_column` or `preliminary_posterior` there.

# Column layout
 - `base_residual_variance`, `G_pi`, `G_variance`, `G_lambda2`, …: gene `c`
   in column `c`.
 - `base_fixedeffects`: the `p` coefficients of gene 1, then gene 2, …, so
   coefficient `k` of gene `c` is in column `(c - 1) * p + k`. With the
   default `X` (intercept only), column `c` is the intercept of gene `c`.

# Example
```julia
Rvar = read_output("output/base_residual_variance.tsv.gz")   # iterations × m
B = read_output("output/base_fixedeffects.tsv.gz")           # iterations × (p * m)
B3 = B[:, (3 - 1) * p .+ (1:p)]                              # gene 3's coefficients
```
"""
function read_output(path::AbstractString; delim::Char='\t')::Matrix{Float64}
    vals = Float64[]
    ncols = Ref(0)
    nrows = _scan_rows(path, delim, nothing) do row
        ncols[] = length(row)
        append!(vals, row)
    end
    nrows == 0 && return zeros(Float64, 0, 0)
    return permutedims(reshape(vals, ncols[], nrows))
end

"""
    preliminary_column(path, i; delim='\\t')

Extract column `i` of an output file, one value per complete iteration, from
an analysis that may still be running or was stopped. This is the trace of
one parameter, e.g. for trace plots or convergence diagnostics.

`i` can be an integer, which returns a `Vector{Float64}`, or a vector of
integers, which returns a `Matrix{Float64}` with one row per iteration and
one column per entry of `i` (in the given order, duplicates allowed). Give
all wanted columns in one call: each call is one pass over the file.

Incomplete iterations are handled as in `preliminary_posterior`: a last row
without an end-of-line character, an incomplete binary row, a binary header
not yet on disk, and a gzip stream without a trailer are all accepted, and
only complete iterations are returned. Only the requested fields are parsed
in text files, and only the requested values are read into memory from
binary files, so a single trace can be extracted from very wide files.
Values in unrequested columns of text files are not checked for being
numeric.

# Column layout
Each row of an effects file holds the effects of gene 1, then gene 2, and
so on. The effect of predictor `q` on gene `c` is in column
`(c - 1) * p + q`, where `p` is the number of genes for `G_effects` and the
number of marker loci for `Q_effects`. Per-gene files (`G_pi`,
`base_residual_variance`, …) have gene `c` in column `c`.

# Example
```julia
# Trace of the residual variance of gene 3
rv3 = preliminary_column("output/base_residual_variance.tsv.gz", 3)

# Traces of the effects of gene 2 and gene 5 on gene 7, with m genes
tr = preliminary_column("output/G_effects.bin", [(7 - 1) * m + 2, (7 - 1) * m + 5])
size(tr)   # (iterations, 2)
```
"""
function preliminary_column(path::AbstractString, i::Integer; delim::Char='\t')::Vector{Float64}
    return vec(preliminary_column(path, [i]; delim=delim))
end

function preliminary_column(path::AbstractString, i::AbstractVector{<:Integer};
    delim::Char='\t')::Matrix{Float64}
    all(>=(1), i) || throw(ArgumentError("Column indices must be positive, got $i"))

    # Read the sorted, unique columns in one pass, then map back to the given order
    cols = sort!(unique(Int.(i)))
    vals = Float64[]
    nrows = _scan_rows(path, delim, cols) do row
        append!(vals, row)
    end
    traces = permutedims(reshape(vals, length(cols), nrows))
    return traces[:, [searchsortedfirst(cols, j) for j in i]]
end

# --- Public post-analysis functions (auto-detect binary vs TSV) ---

function accuracy(T::AbstractString, E::AbstractString)
        Tmat = readdlm(T, ',', Float64, '\n')
        Emat = reshape(posterior_mean(E), size(Tmat))

        return cor(vec(Tmat),vec(Emat))
end

function posterior_mean(path::AbstractString; delim = '\t')
    _is_binary(path) && return _posterior_mean_bin(path)
    
    # Initializing objects
    nlines = 0
    o = Union{Float64,Missing}[]
    
    # Open as gzip if the file ends in .gz, else plain text
    io = _open_output(path)
    
    try
        # Scrolling the file
        for line in eachline(io)
            isempty(strip(line)) && continue # Skip the line if it is empty
            values = tryparse.(Float64, split(strip(line), delim))
            any(isnothing, values) && error("Non-numeric value on line $(nlines + 1)")
            if nlines == 0
                o = Float64.(values)
            else
                o .+= values
            end
            nlines += 1
        end
    finally
        close(io)
    end
    nlines == 0 && error("No data rows found in $path")
    return o ./ nlines
end

"""
    posterior_inclusion_probability(path; delim='\\t')

Compute the exact posterior inclusion probability (PIP) for each effect
from an effects file written by `grn` with a spike-and-slab prior. In each
MCMC iteration, an excluded effect is written as exactly `0.0` and an
included effect is non-zero, so the PIP is the proportion of rows where
the value at each position is non-zero.

This gives the exact PIP — identical to averaging the indicator file —
since indicator files are not written by `grn`.

# Example
```julia
# After running grn:
pip = GRN.posterior_inclusion_probability("output/G_effects.tsv.gz")
# Reshape to gene x gene matrix:
G_pip = reshape(pip, n_genes, n_genes)
# Effects included in >50% of iterations:
G_strong = G_pip .> 0.5
```
"""
function posterior_inclusion_probability(path::AbstractString; delim::Char='\t')
    _is_binary(path) && return _pip_bin(path)

    nlines = 0
    counts = Int[]

    io = _open_output(path)

    try
        for line in eachline(io)
            isempty(strip(line)) && continue
            values = tryparse.(Float64, split(strip(line), delim))
            any(isnothing, values) && error("Non-numeric value on line $(nlines + 1)")
            vals = Float64.(values)
            if nlines == 0
                counts = zeros(Int, length(vals))
            end
            @inbounds for i in eachindex(vals)
                vals[i] != 0.0 && (counts[i] += 1)
            end
            nlines += 1
        end
    finally
        close(io)
    end
    nlines == 0 && error("No data rows found in $path")
    return counts ./ nlines
end

"""
    derive_indicator(effects_path; threshold=0.5, delim='\\t')

Derive a binary indicator vector from an effects file written by `grn`
with a spike-and-slab prior. An indicator is `true` wherever the posterior
inclusion probability (PIP) exceeds `threshold` (default: 0.5, i.e. the
effect was included in more than half of the MCMC iterations).

This uses the exact PIP computed from non-zero counts in the effects file
(`posterior_inclusion_probability`), so it gives the same result as
thresholding the averaged indicator file — without the indicator file
being written.

# Example
```julia
# Effects included in >50% of iterations:
ind = GRN.derive_indicator("output/G_effects.tsv.gz")
# Effects included in >80% of iterations:
ind_strong = GRN.derive_indicator("output/G_effects.tsv.gz"; threshold=0.8)
```
"""
function derive_indicator(effects_path::AbstractString; threshold::Float64=0.5, delim::Char='\t')
    pip = posterior_inclusion_probability(effects_path; delim=delim)
    return pip .> threshold
end
