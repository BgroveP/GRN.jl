# Effect priors for gene-gene (G) and marker-gene (Q) effects
#
# A prior consists of
#  1. A selector type that the user constructs (e.g. `SpikeSlab()`),
#     carrying default settings for the hyperprior
#  2. A state type holding the per-gene Gibbs state, created with `_initstate`
#  3. Methods for `_sample_effects_gene!`, `_sample_hyperparameters!`,
#     `_outputs`, and `_varscale!`
#  4. For continuing an analysis (`grn(path, ...)`): a `_restore!` method if part of the sampled state
#     is not written to file, and a case in `_detect_prior`
#
# Adding a new prior: define a selector and a state subtype of
# `AbstractEffectPriorState`, plus the methods above.

abstract type AbstractEffectPrior end
abstract type AbstractEffectPriorState end

# ---------------------------------------------------------------- Selectors

"""

    SpikeSlab(; includeprob, var, vardf)

Spike-and-slab prior for gene-gene or marker-gene effects. Each effect is
either zero (spike) with probability `1 - pi` or drawn from a normal
distribution with variance `var` (slab). `pi` and `var` are sampled per gene
with `vardf` degrees of freedom for the scaled inverse chi-squared prior on
`var`. Construct as a value in the `priors` NamedTuple, e.g.
`priors = (Gprior = SpikeSlab(),)`.

"""
struct SpikeSlab <: AbstractEffectPrior
    includeprob::Float64
    var::Float64
    vardf::Float64
end
SpikeSlab(; includeprob::Float64=0.5, var::Float64=1.0, vardf::Float64=4.0) =
    SpikeSlab(includeprob, var, vardf)

"""

    Ridge(; var, vardf)

Gaussian prior for all effects of a gene (ridge-like shrinkage, no spike).
The prior variance `var` is sampled per gene with `vardf` degrees of freedom
for the scaled inverse chi-squared prior. Construct as a value in the
`priors` NamedTuple, e.g. `priors = (Gprior = Ridge(),)`.

"""
struct Ridge <: AbstractEffectPrior
    var::Float64
    vardf::Float64
end
Ridge(; var::Float64=1.0, vardf::Float64=4.0) = Ridge(var, vardf)

"""

    Laplace(; lambda2, a, b)

Bayesian lasso prior (Laplace) via a scale mixture of normals: each effect
has its own local variance `tau2` with an exponential prior with rate
`lambda2 / 2`, and `lambda2` per gene has a `Gamma(a, b)` (shape, rate)
hyperprior. Construct as a value in the `priors` NamedTuple, e.g.
`priors = (Gprior = Laplace(),)`.

"""
struct Laplace <: AbstractEffectPrior
    lambda2::Float64
    a::Float64
    b::Float64
end
Laplace(; lambda2::Float64=1.0, a::Float64=1.0, b::Float64=1e-3) =
    Laplace(lambda2, a, b)

# ---------------------------------------------------------------- States

struct SpikeSlabState <: AbstractEffectPriorState
    effects::Vector{Vector{Float64}}
    indicator::Vector{Vector{Int8}}
    pi::Vector{Float64}
    var::Vector{Float64}
    vardf::Vector{Float64}
    varscale::Vector{Float64}
    candidates::Vector{Vector{Int}}
end

struct RidgeState <: AbstractEffectPriorState
    effects::Vector{Vector{Float64}}
    var::Vector{Float64}
    vardf::Vector{Float64}
    varscale::Vector{Float64}
    candidates::Vector{Vector{Int}}
end

struct LaplaceState <: AbstractEffectPriorState
    effects::Vector{Vector{Float64}}
    tau2::Vector{Vector{Float64}}
    lambda2::Vector{Float64}
    a::Vector{Float64}
    b::Vector{Float64}
    candidates::Vector{Vector{Int}}
end

# ---------------------------------------------------------------- Init

function _candidates(genes::Int, predictors::Int; diag::Bool=true)::Vector{Vector{Int}}
    o = [collect(1:predictors) for _ in 1:genes]
    if !diag
        for c in 1:genes
            deleteat!(o[c], c)
        end
    end
    return o
end

function _initstate(sel::SpikeSlab, genes::Int, predictors::Int; diag::Bool=true)
    return SpikeSlabState(
        [zeros(predictors) for _ in 1:genes],
        [zeros(Int8, predictors) for _ in 1:genes],
        fill(sel.includeprob, genes),
        fill(sel.var, genes),
        fill(sel.vardf, genes),
        fill(sel.var * (sel.vardf - 2.0), genes),
        _candidates(genes, predictors; diag=diag))
end

function _initstate(sel::Ridge, genes::Int, predictors::Int; diag::Bool=true)
    return RidgeState(
        [zeros(predictors) for _ in 1:genes],
        fill(sel.var, genes),
        fill(sel.vardf, genes),
        fill(sel.var * (sel.vardf - 2.0), genes),
        _candidates(genes, predictors; diag=diag))
end

function _initstate(sel::Laplace, genes::Int, predictors::Int; diag::Bool=true)
    return LaplaceState(
        [zeros(predictors) for _ in 1:genes],
        [fill(1.0 / sel.lambda2, predictors) for _ in 1:genes],
        fill(sel.lambda2, genes),
        fill(sel.a, genes),
        fill(sel.b, genes),
        _candidates(genes, predictors; diag=diag))
end

# ---------------------------------------------------------------- Outputs

_outputs(st::SpikeSlabState) = (
    effects=st.effects,
    pi=st.pi,
    variance=st.var)

_outputs(st::RidgeState) = (
    effects=st.effects,
    variance=st.var)

_outputs(st::LaplaceState) = (
    effects=st.effects,
    tau2=st.tau2,
    lambda2=st.lambda2)

# Prior variance scale for the scaled inverse chi-squared variance priors
_varscale!(st::SpikeSlabState) = st.varscale .= st.var .* (st.vardf .- 2.0)
_varscale!(st::RidgeState) = st.varscale .= st.var .* (st.vardf .- 2.0)
_varscale!(st::LaplaceState) = nothing

# ---------------------------------------------------------------- Continuation

"""

    _restore!(st)

Rebuild the parts of a prior state that are sampled but not written to the
output files, after the outputs (`_outputs(st)`) have been read back from
the last stored iteration. Used when `grn` continues an analysis. For spike-and-slab, the
indicators are recovered from the effects (an included effect is non-zero
with probability one). Priors whose state is fully written need no method.

"""
_restore!(::AbstractEffectPriorState)::Nothing = nothing
function _restore!(st::SpikeSlabState)::Nothing
    for c in eachindex(st.effects)
        st.indicator[c] .= st.effects[c] .!= 0.0
    end
    return nothing
end

"""

    _detect_prior(path, group, ext)

Return the type of the prior selector used for `group` (`:G` or `:Q`) in an
earlier analysis, from the output files in `path` with extension `ext`:
`<group>_pi` exists for `SpikeSlab`, `<group>_tau2` for `Laplace`, and
`<group>_variance` alone for `Ridge`. Returns `nothing` when the group has
no output files.

"""
function _detect_prior(path::AbstractString, group::Symbol,
    ext::AbstractString)::Union{Nothing, Type{<:AbstractEffectPrior}}
    has(name) = isfile(joinpath(path, "$(group)_$name$ext"))
    has("effects") || return nothing
    has("pi") && return SpikeSlab
    has("tau2") && return Laplace
    has("variance") && return Ridge
    return nothing
end
