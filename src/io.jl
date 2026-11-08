# Functions for creating output files and writing parameters to them

# Binary format magic bytes ("GRNB") and version
const _BIN_MAGIC = UInt32(0x42524E47)  # "GRNB" in little-endian
const _BIN_VERSION = UInt32(1)
const _BIN_HEADER_SIZE = 16  # magic(4) + version(4) + ncols(8)

"""
    makeout(p, outputs; force=false, compress=true, format=:tsv, append=false)

Create an output folder and open one output file per entry in
`outputs`, a NamedTuple that maps a file name to the object that will be
written to it. Returns a vector of NamedTuples with the name, the output
object, the IO stream, a persistent IOBuffer, and the format.

# Arguments
 - `p`:        Output directory path
 - `outputs`:  NamedTuple mapping file names to output objects
 - `force`:    Overwrite existing directory (ignored with `append=true`)
 - `compress`: Gzip-compress the output (level 1, fast)
 - `format`:   `:tsv` for tab-separated text (default) or `:bin` for
               raw binary Float64 with a 16-byte header
 - `append`:   Open the existing files in `p` and append rows to them
               instead of creating new files. Every file must exist. No
               binary header is written, and gzip output is added as a new
               gzip member, which gzip readers decompress as one stream.
               Used when `grn` continues an analysis.

# File extensions
 - `format=:tsv, compress=true`  → `.tsv.gz`
 - `format=:tsv, compress=false` → `.tsv`
 - `format=:bin, compress=true`  → `.bin.gz`
 - `format=:bin, compress=false` → `.bin`
"""
function makeout(p::AbstractString, outputs::NamedTuple; force::Bool=false,
    compress::Bool=true, format::Symbol=:tsv, append::Bool=false)::Vector{NamedTuple}

    format in (:tsv, :bin) || throw(ArgumentError(
        "format must be :tsv or :bin, got :$format"))

    if append
        # Append to the files of an existing directory
        isdir(p) || error("The output directory $p does not exist")
    else
        # Handle existing directory
        if isdir(p)
            force || error("The output directory $p already exists")
            rm(p, force=true, recursive=true)
            sleep(1) # Wait for the removal to be picked up by the file system
        end

        # Make new directory
        mkpath(p)
    end

    # Determine file extension
    ext = if format === :bin
        compress ? ".bin.gz" : ".bin"
    else
        compress ? ".tsv.gz" : ".tsv"
    end

    # Check that all files exist before any of them is opened for appending
    if append
        for k in keys(outputs)
            f = joinpath(abspath(p), "$k$ext")
            isfile(f) || error("Cannot append: the output file $f does not exist")
        end
    end

    # Create one IO stream and one persistent buffer per output object
    mode = append ? "a" : "w"
    files = NamedTuple[]
    for k in keys(outputs)
        f = joinpath(abspath(p), "$k$ext")
        io = compress ? GZip.open(f, mode, 1) : open(f, mode)
        buf = IOBuffer(sizehint=_BUFSIZE)

        # Write binary header (magic + version + ncols) to new files only
        if format === :bin && !append
            ncols = _ncols(outputs[k])
            write(buf, _BIN_MAGIC)
            write(buf, _BIN_VERSION)
            write(buf, Int64(ncols))
        end

        push!(files, (name=k, values=outputs[k], io=io, buf=buf, format=format))
    end
    return files
end

"""
    _ncols(x)

Return the number of Float64 values in one row of output object `x`.
"""
_ncols(x::AbstractVector{Float64})::Int = length(x)
_ncols(x::AbstractVector{<:AbstractVector{Float64}})::Int = sum(length, x)
_ncols(x::AbstractVector{<:Real})::Int = length(x)
_ncols(x::AbstractVector{<:AbstractVector})::Int = sum(length, x)

"""
    write_output(o; sep='\t')

Write one row of the output object in `o` to its IOBuffer. For `format=:tsv`,
vectors are written as one row of tab-separated values, and vectors of vectors
are written as one flattened row. For `format=:bin`, raw Float64 bytes are
written with no separators. Uses a persistent IOBuffer (1 MB, created in
`makeout`) that is flushed to the stream whenever it fills, so memory usage
stays bounded even for very large rows (e.g. 17K × 17K effects).

The buffer is NOT flushed after every row — data accumulates across rows
and is only written to the stream when the buffer fills. Call `flush_outputs`
before closing the streams to write any remaining buffered data.
"""
function write_output(o::NamedTuple; sep::Char='\t')::Nothing
    if get(o, :format, :tsv) === :bin
        _write_row_bin!(o.buf, o.io, o.values)
    else
        _write_row!(o.buf, o.io, o.values, sep)
    end
    return nothing
end

# Flush the buffer to the stream if it has reached the size threshold.
@inline function _flush_if_full!(buf::IOBuffer, io::IO)::Nothing
    if position(buf) >= _BUFSIZE
        write(io, take!(buf))
    end
    return nothing
end

# --- Fast path: 1D Float64 vector ---

function _write_row!(buf::IOBuffer, io::IO, x::AbstractVector{Float64}, sep::Char)::Nothing
    @inbounds for i in eachindex(x)
        i != firstindex(x) && print(buf, sep)
        @printf(buf, "%.10g", x[i])
        _flush_if_full!(buf, io)
    end
    println(buf)
    return nothing
end

# --- Fast path: nested Float64 vectors (flattened into one row) ---

function _write_row!(buf::IOBuffer, io::IO, x::AbstractVector{<:AbstractVector{Float64}}, sep::Char)::Nothing
    first = true
    for v in x
        @inbounds for val in v
            if first
                @printf(buf, "%.10g", val)
                first = false
            else
                print(buf, sep)
                @printf(buf, "%.10g", val)
            end
        end
        _flush_if_full!(buf, io)
    end
    println(buf)
    return nothing
end

# --- Generic fallbacks for non-Float64 element types ---

function _write_row!(buf::IOBuffer, io::IO, x::AbstractVector{<:Real}, sep::Char)::Nothing
    first = true
    for val in x
        first || print(buf, sep)
        print(buf, val)
        first = false
        _flush_if_full!(buf, io)
    end
    println(buf)
    return nothing
end

function _write_row!(buf::IOBuffer, io::IO, x::AbstractVector{<:AbstractVector}, sep::Char)::Nothing
    first = true
    for v in x
        for val in v
            first || print(buf, sep)
            print(buf, val)
            first = false
            _flush_if_full!(buf, io)
        end
    end
    println(buf)
    return nothing
end

# --- Binary output: raw Float64 bytes, no separators ---
#
# Binary writes use unsafe_write with pointer arithmetic for zero-copy,
# allocation-free streaming. Data is written in chunks of _BUFSIZE bytes
# so the IOBuffer never grows beyond ~1 MB regardless of row size.

"""
    _write_row_bin!(buf, io, x)

Write one row of Float64 values to `buf` in raw binary format (native
endianness, no separators). Flushes to `io` when the buffer fills.
"""
function _write_row_bin!(buf::IOBuffer, io::IO, x::AbstractVector{Float64})::Nothing
    n = length(x)
    n == 0 && return nothing
    ptr = Ptr{UInt8}(pointer(x))
    total_bytes = n * sizeof(Float64)
    offset = 0
    while offset < total_bytes
        chunk = min(_BUFSIZE, total_bytes - offset)
        unsafe_write(buf, ptr + offset, chunk)
        offset += chunk
        _flush_if_full!(buf, io)
    end
    return nothing
end

function _write_row_bin!(buf::IOBuffer, io::IO, x::AbstractVector{<:AbstractVector{Float64}})::Nothing
    for v in x
        _write_row_bin!(buf, io, v)
    end
    return nothing
end

# Generic fallbacks for non-Float64 element types (convert to Float64)
function _write_row_bin!(buf::IOBuffer, io::IO, x::AbstractVector{<:Real})::Nothing
    _write_row_bin!(buf, io, Float64.(x))
    return nothing
end

function _write_row_bin!(buf::IOBuffer, io::IO, x::AbstractVector{<:AbstractVector})::Nothing
    for v in x
        _write_row_bin!(buf, io, Float64.(v))
    end
    return nothing
end

function _print(i::NamedTuple)::Nothing
    println("Number of individuals:   $(i.individuals)")
    println("Number of genes:         $(i.genes)")
    println("Number of markers:       $(i.loci)")
    println("Number of fixed effects: $(get(i, :fixed, 1))")
    println("Number of threads:       $(nthreads())")
    return nothing
end

"""
    flush_outputs(outs)

Flush all pending IOBuffer data to the underlying streams. Must be called
before `close_outputs` (or before closing streams manually) to ensure no
buffered data is lost.
"""
function flush_outputs(outs::Vector{NamedTuple})::Nothing
    for o in outs
        write(o.io, take!(o.buf))
    end
    return nothing
end

"""
    close_outputs(outs)

Flush all IOBuffer data and close all IO streams in the vector of output
NamedTuples. Must be called to flush gzip trailers and release file handles.
"""
function close_outputs(outs::Vector{NamedTuple})::Nothing
    flush_outputs(outs)
    for o in outs
        close(o.io)
    end
    return nothing
end
