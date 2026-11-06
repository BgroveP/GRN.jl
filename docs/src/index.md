# GRN

`GRN.jl` estimates gene regulatory networks from gene expression data,
optionally jointly with genotype (marker) data, using a Bayesian linear
model and Gibbs sampling.

- **Repository**: <https://github.com/BgroveP/GRN>
- **Getting started**: see the README for installation and a first run.

## Priors

The prior for each effect group is selectable per relation type:

- `SpikeSlab(; includeprob, var, vardf)` — spike-and-slab prior (default)
- `Ridge(; var, vardf)` — Gaussian prior with sampled per-gene variance
- `Laplace(; lambda2, a, b)` — Bayesian lasso prior

## Post-analysis

Use the functions in `GRN.postanalysis` (`accuracy`, `posterior_mean`,
`posterior_inclusion_probability`, `derive_indicator`) to summarize the
sampled output files.

`preliminary_posterior(path)` returns the posterior mean, PIP, and number of
complete iterations from an output file of an analysis that is still running
or was stopped, ignoring a last row that was cut off mid-write.
