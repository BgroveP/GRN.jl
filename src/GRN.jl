module GRN

# Used packages
using Base.Threads
using LinearAlgebra
using Printf
using Random
using Statistics
using DelimitedFiles
using GZip
using Zlib_jll: libz

# Set number of BLAS threads to one
BLAS.set_num_threads(1)

# Fixed buffer size for streaming output writes (1 MB).
# Bounds memory per output file regardless of row size.
const _BUFSIZE = 1 << 20

# Includes
include("postanalysis.jl")
include("priors.jl")
include("initialize.jl")
include("check.jl")
include("estimate.jl")
include("io.jl")
include("sample.jl")
include("simulate.jl")
include("utility.jl")

# Exports
export grn, simulate
export SpikeSlab, Ridge, Laplace
export accuracy, posterior_mean, posterior_inclusion_probability, derive_indicator
export read_binary, preliminary_posterior

end
