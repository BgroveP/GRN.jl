# MCMC convergence diagnostics
#
# Effective sample size following Vehtari, Gelman, Simpson, Carpenter and
# Bürkner (2021), "Rank-normalization, folding, and localization: An
# improved R-hat for assessing convergence of MCMC", Bayesian Analysis
# 16(2), 667-718 — the method used by Stan and the R package `posterior`.
#
# Adapted from posterior (BSD-3-Clause), https://github.com/stan-dev/posterior,
# effective sample size functions (ess_bulk, ess_tail, ess_basic, .ess).
# Copyright and license text: THIRD_PARTY_NOTICES.md.
#
# Chains are given as a matrix with one column per chain (iterations ×
# chains). A single trace is one chain.

"""
    _qnorm(p)

Quantile function of the standard normal distribution (Wichura 1988,
algorithm AS 241, accurate to about 1e-16). Used for rank normalization, so
`Distributions` is not needed.
"""
function _qnorm(p::Float64)::Float64
    0.0 < p < 1.0 || return p == 0.0 ? -Inf : p == 1.0 ? Inf : NaN
    q = p - 0.5
    if abs(q) <= 0.425
        r = 0.180625 - q * q
        return q * (((((((r * 2509.0809287301226727 +
                          33430.575583588128105) * r + 67265.770927008700853) * r +
                        45921.953931549871457) * r + 13731.693765509461125) * r +
                      1971.5909503065514427) * r + 133.14166789178437745) * r +
                    3.387132872796366608) /
               (((((((r * 5226.495278852545925 +
                      28729.085735721942674) * r + 39307.89580009271061) * r +
                    21213.794301586595867) * r + 5394.1960214247511077) * r +
                  687.1870074920579083) * r + 42.313330701600911252) * r + 1.0)
    end
    r = sqrt(-log(q < 0.0 ? p : 1.0 - p))
    if r <= 5.0
        r -= 1.6
        val = (((((((r * 7.7454501427834140764e-4 +
                     0.0227238449892691845833) * r + 0.24178072517745061177) * r +
                   1.27045825245236838258) * r + 3.64784832476320460504) * r +
                 5.7694972214606914055) * r + 4.6303378461565452959) * r +
               1.42343711074968357734) /
              (((((((r * 1.05075007164441684324e-9 +
                     5.475938084995344946e-4) * r + 0.0151986665636164571966) * r +
                   0.14810397642748007459) * r + 0.68976733498510000455) * r +
                 1.6763848301838038494) * r + 2.05319162663775882187) * r + 1.0)
    else
        r -= 5.0
        val = (((((((r * 2.01033439929228813265e-7 +
                     2.71155556874348757815e-5) * r + 0.0012426609473880784386) * r +
                   0.026532189526576123093) * r + 0.29656057182850489123) * r +
                 1.7848265399172913358) * r + 5.4637849111641143699) * r +
               6.6579046435011037772) /
              (((((((r * 2.04426310338993978564e-15 +
                     1.4215117583164458887e-7) * r + 1.8463183175100546818e-5) * r +
                   7.868691311456132591e-4) * r + 0.0148753612908506148525) * r +
                 0.13692988092273580531) * r + 0.59983220655588793769) * r + 1.0)
    end
    return q < 0.0 ? -val : val
end

"""
    _split_chains(x)

Split each chain (column) into its first and second half, giving twice as
many chains of half the length. For an odd number of iterations the middle
iteration is dropped.
"""
function _split_chains(x::AbstractMatrix{Float64})::Matrix{Float64}
    n = size(x, 1)
    half = n ÷ 2
    return hcat(x[1:half, :], x[n-half+1:n, :])
end

"""
    _rank_normalize(x)

Replace all values of `x` by their normal scores: the average rank `r` over
all chains (ties share their average rank) mapped through
`qnorm((r - 3/8) / (S + 1/4))`, with `S` the total number of draws.
"""
function _rank_normalize(x::AbstractMatrix{Float64})::Matrix{Float64}
    S = length(x)
    p = sortperm(vec(x))
    z = similar(x)
    i = 1
    while i <= S
        # Find the run of tied values starting at position i
        j = i
        while j < S && x[p[j+1]] == x[p[i]]
            j += 1
        end
        zr = _qnorm(((i + j) / 2 - 0.375) / (S + 0.25))
        for k in i:j
            z[p[k]] = zr
        end
        i = j + 1
    end
    return z
end

"""
    _autocov(x, means, t)

Autocovariance at lag `t` averaged over chains (columns of `x`), with the
biased estimator `sum((x[i] - m) * (x[i+t] - m)) / n` per chain.
"""
function _autocov(x::AbstractMatrix{Float64}, means::Vector{Float64}, t::Int)::Float64
    n, m = size(x)
    s = 0.0
    @inbounds for j in 1:m
        mu = means[j]
        cj = 0.0
        for i in 1:n-t
            cj += (x[i, j] - mu) * (x[i+t, j] - mu)
        end
        s += cj / n
    end
    return s / m
end

"""
    _ess(x)

Effective sample size of the draws in `x` (iterations × chains), already
split and transformed as wanted. Autocorrelations are combined over chains
and truncated with Geyer's initial positive sequence, then made monotone
with Geyer's initial monotone sequence. Returns `NaN` for fewer than 5
iterations per chain, non-finite values, or constant draws.
Autocovariances are computed only up to the truncation lag, so well-mixing
chains are cheap.
"""
function _ess(x::AbstractMatrix{Float64})::Float64
    n, m = size(x)
    (n < 5 || any(!isfinite, x)) && return NaN
    x1 = x[1]
    all(==(x1), x) && return NaN

    # Within-chain variance and the combined variance estimate
    means = vec(mean(x, dims=1))
    mean_var = _autocov(x, means, 0) * n / (n - 1)
    var_plus = mean_var * (n - 1) / n
    m > 1 && (var_plus += var(means))

    # rho[t + 1] is the autocorrelation at lag t
    rho = zeros(Float64, n)
    rho[1] = 1.0
    rho_even = 1.0
    rho_odd = 1.0 - (mean_var - _autocov(x, means, 1)) / var_plus
    rho[2] = rho_odd

    # Geyer's initial positive sequence: stop when a pair of lags sums to < 0
    t = 0
    while t < n - 5 && !isnan(rho_even + rho_odd) && rho_even + rho_odd > 0.0
        t += 2
        rho_even = 1.0 - (mean_var - _autocov(x, means, t)) / var_plus
        rho_odd = 1.0 - (mean_var - _autocov(x, means, t + 1)) / var_plus
        if rho_even + rho_odd >= 0.0
            rho[t+1] = rho_even
            rho[t+2] = rho_odd
        end
    end
    max_t = t
    rho_even > 0.0 && (rho[max_t+1] = rho_even)

    # Geyer's initial monotone sequence: pair sums must not increase
    t = 0
    while t <= max_t - 4
        t += 2
        if rho[t+1] + rho[t+2] > rho[t-1] + rho[t]
            rho[t+1] = (rho[t-1] + rho[t]) / 2.0
            rho[t+2] = rho[t+1]
        end
    end

    # Integrated autocorrelation time, bounded so ESS <= S * log10(S)
    S = n * m
    tau = -1.0 + 2.0 * sum(@view rho[1:max_t]) + rho[max_t+1]
    tau = max(tau, 1.0 / log10(S))
    return S / tau
end

"""
    ess(x; kind=:bulk)

Effective sample size of MCMC draws: the number of independent draws that
would estimate the posterior mean as precisely as the correlated draws in
`x`. Gibbs draws are autocorrelated, so the ESS is usually well below the
number of iterations; posterior means are reliable when the ESS is at least
a few hundred.

`x` is a vector (one chain) or a matrix with one column per chain
(iterations × chains, equal lengths). The chains are split in half before
the computation, so a chain that drifts gives a low ESS.

# Kinds
 - `:bulk`:  ESS of the rank-normalized draws (default). Robust to heavy
             tails and to the many exact zeros of spike-and-slab effects;
             describes the centre of the distribution (mean, median).
 - `:tail`:  Minimum of the ESS of the indicators `x <= q05` and `x <= q95`,
             with `q05` and `q95` the 5% and 95% quantiles of all draws;
             describes the reliability of the 90% interval.
 - `:basic`: ESS of the draws without rank normalization.

Returns `NaN` when the ESS is undefined: fewer than 10 iterations per chain
(5 per split half), non-finite values, or constant draws (e.g. an effect
that is zero in every iteration). For `:tail`, `NaN` is also returned when
an indicator is constant.

Method: Vehtari et al. (2021), Bayesian Analysis 16(2), 667-718, as in Stan
and the R package `posterior`. Results agree with
`MCMCDiagnosticTools.ess(...; maxlag=typemax(Int))` to rounding error; that
package by default truncates the autocorrelations at lag 250, which gives
higher ESS for very slowly mixing chains.

# Example
```julia
tr = preliminary_column("output/base_residual_variance.tsv.gz", 3)
ess(tr)                  # bulk ESS
ess(tr; kind=:tail)

# ESS of the inclusion indicator, which the PIP is the mean of
w = preliminary_column("output/G_effects.tsv.gz", (7 - 1) * m + 2)
ess(w .!= 0; kind=:basic)
```
"""
function ess(x::AbstractMatrix{<:Real}; kind::Symbol=:bulk)::Float64
    kind in (:bulk, :tail, :basic) || throw(ArgumentError(
        "kind must be :bulk, :tail, or :basic, got :$kind"))
    xf = Matrix{Float64}(x)
    (isempty(xf) || any(!isfinite, xf)) && return NaN

    # Ranks and quantiles use all draws, including a middle iteration that
    # the split drops for an odd number of iterations (as MCMCDiagnosticTools)
    if kind === :bulk
        return _ess(_split_chains(_rank_normalize(xf)))
    elseif kind === :basic
        return _ess(_split_chains(xf))
    else
        xs = _split_chains(xf)
        q05, q95 = quantile(vec(xf), (0.05, 0.95))
        return min(_ess(Float64.(xs .<= q05)), _ess(Float64.(xs .<= q95)))
    end
end

ess(x::AbstractVector{<:Real}; kind::Symbol=:bulk)::Float64 =
    ess(reshape(x, :, 1); kind=kind)

"""
    ess(paths, i; kind=:bulk, delim='\\t')

Effective sample size of column `i` of one or more output files, read with
`preliminary_column`, so files of a running analysis can be used. `paths`
is one file path, or a vector of paths to the same output file from
independent runs (one chain each, e.g. different seeds and `outpath`s).
Chains are truncated to the shortest one.

`i` is an integer (returns a `Float64`) or a vector of integers (returns a
`Vector{Float64}`, one ESS per column, from one pass over each file). See
`preliminary_column` for the column layout of effects files.

# Example
```julia
ess("output/G_pi.tsv.gz", 1:m)                           # every gene's pi
ess(["run1/G_pi.bin", "run2/G_pi.bin", "run3/G_pi.bin"], 4)  # three chains
```
"""
function ess(paths::AbstractVector{<:AbstractString}, i::AbstractVector{<:Integer};
    kind::Symbol=:bulk, delim::Char='\t')::Vector{Float64}
    isempty(paths) && throw(ArgumentError("No files given"))
    traces = [preliminary_column(p, i; delim=delim) for p in paths]
    n = minimum(size(t, 1) for t in traces)
    return [ess(reduce(hcat, (t[1:n, k] for t in traces)); kind=kind) for k in eachindex(i)]
end

ess(paths::AbstractVector{<:AbstractString}, i::Integer; kind::Symbol=:bulk,
    delim::Char='\t')::Float64 = only(ess(paths, [i]; kind=kind, delim=delim))

ess(path::AbstractString, i::Union{Integer, AbstractVector{<:Integer}};
    kind::Symbol=:bulk, delim::Char='\t') = ess([path], i; kind=kind, delim=delim)
