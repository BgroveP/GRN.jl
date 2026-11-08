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

## Continuing an analysis

`grn(path, Y, M, chain)` continues an analysis from the last
iteration stored in the output folder `path`, and appends the new iterations
to the existing files. Give the same `Y`, `M`, `X`, and `priors` as in the
call to `grn`; the prior types, format, and compression are read from the
files. An analysis that was interrupted (unfinished files) is continued from
the last iteration complete in all files, into a new folder given with
`outpath`.

## Post-analysis

Use the functions in `GRN.postanalysis` (`accuracy`, `posterior_mean`,
`posterior_inclusion_probability`, `derive_indicator`) to summarize the
sampled output files.

`preliminary_posterior(path)` returns the posterior mean, PIP, and number of
complete iterations from an output file of an analysis that is still running
or was stopped, ignoring a last row that was cut off mid-write.
`preliminary_column(path, i)` extracts the trace of column `i` (or several
columns) from such files, and `read_output(path)` reads all iterations of a
small output file (e.g. residual variances, fixed effects) into a matrix.

## Convergence diagnostics

`ess(x; kind=:bulk)` computes the effective sample size of a trace or of
several chains (Vehtari et al. 2021), and `ess(paths, i)` does so directly
from output files.
