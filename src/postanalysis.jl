
        
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
    _read_partial!(src, buf)

Read bytes into `buf` until it is full or the end of the available data is
reached. Returns the number of bytes read; fewer than `length(buf)` means the
end of the data written so far. A gzip stream that ends without a trailer is
treated as the end of the data, not as an error.
"""
function _read_partial!(src::IOStream, buf::Vector{UInt8})::Int
    return readbytes!(src, buf, length(buf))
end

function _read_partial!(src::_PartialGz, buf::Vector{UInt8})::Int
    nr = 0
    while nr < length(buf)
        # gzread returns a C int, so read at most 1 GiB per call
        len = min(length(buf) - nr, 1 << 30)
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
    _preliminary_bin(path)

Posterior summaries from the complete rows of a binary output file that may
still be written. A missing or incomplete header gives zero rows, and an
incomplete last row is ignored.
"""
function _preliminary_bin(path::AbstractString)::NamedTuple
    src = _open_partial(path)
    try
        hdr = Vector{UInt8}(undef, _BIN_HEADER_SIZE)
        _read_partial!(src, hdr) < _BIN_HEADER_SIZE &&
            return _preliminary_result(Float64[], Int[], 0)
        ncols = _read_bin_header(IOBuffer(hdr))

        rowbytes = Vector{UInt8}(undef, ncols * sizeof(Float64))
        row = reinterpret(Float64, rowbytes)
        sums = zeros(Float64, ncols)
        counts = zeros(Int, ncols)
        nrows = 0
        while _read_partial!(src, rowbytes) == length(rowbytes)
            _accumulate!(sums, counts, row)
            nrows += 1
        end
        return _preliminary_result(sums, counts, nrows)
    finally
        _close_partial(src)
    end
end

"""
    _preliminary_tsv(path, delim)

Posterior summaries from the complete rows of a text output file that may
still be written. Only rows terminated by an end-of-line character are used;
an unterminated last row (possibly cut mid-number) is ignored. The file is
parsed in 1 MB chunks, so memory use is one row of values plus the
accumulators, not one row of text. Blank lines and `\\r\\n` line endings are
accepted.
"""
function _preliminary_tsv(path::AbstractString, delim::Char)::NamedTuple
    isascii(delim) && delim != '\n' || throw(ArgumentError(
        "delim must be an ASCII character other than newline, got $(repr(delim))"))
    d = UInt8(delim)
    nl = UInt8('\n')

    src = _open_partial(path)
    chunk = Vector{UInt8}(undef, _BUFSIZE)
    carry = UInt8[]     # bytes of a field that continues in the next chunk
    row = Float64[]     # values of the current row
    sums = Float64[]
    counts = Int[]
    ncols = 0           # values per row, set by the first complete row
    col = 0             # values read so far in the current row
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

                v = tryparse(Float64, field)
                if v === nothing
                    # Skip blank lines; anything else is malformed
                    b == nl && col == 0 && all(isspace, field) && continue
                    error("Non-numeric value $(repr(field)) in row $(nrows + 1) of $path")
                end

                col += 1
                if nrows == 0
                    push!(row, v)
                else
                    col <= ncols || error(
                        "Row $(nrows + 1) of $path has more than $ncols values")
                    @inbounds row[col] = v
                end

                # End of line: the row is complete
                if b == nl
                    if nrows == 0
                        ncols = col
                        sums = zeros(Float64, ncols)
                        counts = zeros(Int, ncols)
                    else
                        col == ncols || error(
                            "Row $(nrows + 1) of $path has $col values, expected $ncols")
                    end
                    _accumulate!(sums, counts, row)
                    nrows += 1
                    col = 0
                end
            end
            append!(carry, view(chunk, start:n))
        end
    finally
        _close_partial(src)
    end

    # Values in `row` and bytes in `carry` belong to an unterminated last row
    # and are dropped
    return _preliminary_result(sums, counts, nrows)
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
    return _is_binary(path) ? _preliminary_bin(path) : _preliminary_tsv(path, delim)
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
