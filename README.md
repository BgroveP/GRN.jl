# GRN

[![Build Status](https://github.com/BgroveP/GRN/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/BgroveP/GRN/actions/workflows/CI.yml?query=branch%3Amain)

## Overview

GRN.jl estimates gene regulatory networks from gene expression data, optionally jointly with genotype (marker) data, using a Bayesian linear model and Gibbs sampling. The model links each gene's expression to the expression of the other genes (gene-gene effects) and, when marker data are supplied, to marker genotypes (marker-gene effects).

The prior for each effect group is selectable per relation type:

- `SpikeSlab(; includeprob, var, vardf)` — spike-and-slab prior with sampled inclusion probabilities (default)
- `Ridge(; var, vardf)` — Gaussian prior on all effects with a sampled per-gene prior variance
- `Laplace(; lambda2, a, b)` — Bayesian lasso prior with per-effect local variances and a sampled per-gene penalty.

Sampled parameters are written to tab-separated output files, one row per retained MCMC iteration.

## Requirements

- Julia 1.10 or newer
- The `DelimitedFiles` standard library to load CSV input in the example below

## How to run

Install the package from this repository, for example from a Julia session:

```julia
using Pkg
Pkg.add(url = "https://github.com/BgroveP/GRN")
```

or, for development, clone the repository and activate the package in the project folder with `Pkg.develop(path = "...")`.

The analysis is run with the `grn` function, which takes a matrix `Y` of gene expressions (observations by genes) and, optionally, a matrix `M` of genotypes (observations by markers):

```julia
using DelimitedFiles
using GRN

# Load the input
Y = readdlm("data/Y.csv", ',')   # n × m gene expressions
M = readdlm("data/M.csv", ',')   # n × u genotypes

# Run the MCMC: iterations 1:9999 are burn-in, and parameters are written to file for iterations 10000, 10050, ..., 100000
grn(Y, M, 10000:50:100000, outpath = "grn_example")
```

The `chain` argument (an increasing range of iteration numbers) controls both the number of MCMC iterations run and which of them are written to the output files: iterations before `first(chain)` are treated as burn-in, and every iteration in `chain` is saved. A convenience method `grn(Y, chain; ...)` runs an analysis without genotype data.

### Priors

The priors for gene-gene (`G`) and marker-gene (`Q`) effects are selected with the `Gprior` and `Qprior` keys, and hyperparameters are set with the remaining keys of the `priors` NamedTuple:

```julia
# Spike-and-slab priors with custom hyperparameter priors
grn(Y, M, 10000:50:100000, outpath = "grn_example2", priors = (Gpi = 0.01, Qpi = 0.01, Rvar = 0.01))

# Bayesian lasso prior for gene-gene effects, ridge prior for marker-gene effects
grn(Y, M, 10000:50:100000, outpath = "grn_example3", priors = (Gprior = Laplace(), Qprior = Ridge()))
```

A complete list of prior keys is in the docstring of `grn`.

### Output

The analysis writes one tab-separated file (`.tsv`) per output parameter to the
`outpath` folder, one row per retained iteration. When genotype data are
supplied, files are prefixed with `G_` for gene-gene effects and `Q_` for
marker-gene effects (for example `G_effects.tsv` and `Q_effects.tsv`); without
genotype data only the `G_` and base files are created.

### Preliminary results

`preliminary_posterior` summarizes an output file while the analysis is still
running, or after it was stopped. It uses only complete iterations and
ignores a last row that was cut off mid-write:

```julia
res = preliminary_posterior("grn_example/G_effects.tsv.gz")
res.iterations   # complete iterations on disk so far
res.mean         # posterior mean of each effect
res.pip          # posterior inclusion probability of each effect
```

`preliminary_column` extracts the trace of one or more columns, for example
for trace plots. The effect of predictor `q` on gene `c` is in column
`(c - 1) * p + q`, with `p` the number of genes (`G_effects`) or marker loci
(`Q_effects`):

```julia
rv3 = preliminary_column("grn_example/base_residual_variance.tsv.gz", 3)
tr = preliminary_column("grn_example/G_effects.tsv.gz", [(7 - 1) * m + 2, (7 - 1) * m + 5])
```

Output is buffered, so the files on disk lag somewhat behind the sampler.

### Convergence diagnostics

`ess` computes the effective sample size (bulk, tail, or basic; Vehtari et
al. 2021) from a trace, a matrix of chains, or directly from output files:

```julia
ess("grn_example/G_pi.tsv.gz", 1:m)              # bulk ESS of every gene's pi
ess(["run1/G_pi.bin", "run2/G_pi.bin"], 4)        # two independent chains
ess(preliminary_column("grn_example/G_effects.tsv.gz", 12) .!= 0; kind=:basic)  # inclusion indicator
```

`MCMCDiagnosticTools` also exports `ess`; if both packages are loaded, call
`GRN.ess`.
