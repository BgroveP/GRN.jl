# Changelog

## [Unreleased]

### Added

- `preliminary_posterior(path; delim='\t')` — posterior mean, PIP, and the
  number of complete iterations from an output file of a run that is still
  going or was stopped. Returns `(mean, pip, iterations)` in one pass. Rows
  without an end-of-line character (cut mid-write) and incomplete binary
  rows are ignored; a binary header not yet on disk gives zero iterations;
  gzip streams without a trailer are read up to the last decompressible
  data (GZip.jl throws `ZError(-5)` on these and discards what it read, so
  they are read with zlib's `gzread` directly). Text is parsed in 1 MB
  chunks, so memory is one row of values rather than one row of text. Works
  for `.tsv`, `.tsv.gz`, `.bin`, `.bin.gz`, and `\r\n` line endings
  (`src/postanalysis.jl`)
- `preliminary_column(path, i; delim='\t')` — the trace of column `i`
  (one value per complete iteration) from an output file of a run that is
  still going or was stopped. `i` may be a vector of columns, returning an
  iterations × columns matrix from one pass over the file. Handles unfinished
  files like `preliminary_posterior`. Only the requested fields are parsed in
  text files, and binary files skip unrequested values without holding a
  row in memory. Measured on a 40 × 1M file: 0.09 s (`.bin`), 0.8 s
  (`.bin.gz`), 2.6 s (`.tsv.gz`, vs 9.6 s for `preliminary_posterior`).
  `preliminary_posterior` and `preliminary_column` share new row scanners
  `_scan_tsv_rows` / `_scan_bin_rows` (`src/postanalysis.jl`)
- `read_output(path; delim='\t')` — read all complete iterations of an
  output file into an iterations × values matrix, e.g. for
  `base_residual_variance` and `base_fixedeffects`. Works for all four
  formats and for files of a running or stopped analysis (built on the same
  row scanners as `preliminary_posterior`). Holds the whole file in memory,
  so it is meant for the small per-gene files (`src/postanalysis.jl`)
- `ess(x; kind=:bulk)` — effective sample size of MCMC draws (vector, or
  iterations × chains matrix), following Vehtari et al. (2021): split
  chains, Geyer's initial positive and monotone sequences, and `kind` =
  `:bulk` (rank-normalized, default), `:tail` (5%/95% quantile indicators),
  or `:basic`. File methods `ess(path_or_paths, i)` read column(s) `i` with
  `preliminary_column`, so running analyses and several independent chains
  can be diagnosed directly. Agrees with
  `MCMCDiagnosticTools.ess(...; maxlag=typemax(Int))` to ~1e-14 on the
  tested cases (i.i.d., AR(1), antithetic, ties at zero, heavy tails, odd
  lengths, short chains). No new dependency: the normal quantile function
  for rank normalization is implemented internally (`_qnorm`, AS 241).
  0.16 s for a 100,000-iteration chain with autocorrelation 0.999. The ESS
  computation is adapted from the R package `posterior` (BSD-3-Clause);
  attribution and license text are in `THIRD_PARTY_NOTICES.md`
  (`src/diagnostics.jl`)
- `THIRD_PARTY_NOTICES.md` for attribution and license texts of adapted
  third-party code, and a "Third-party code and licenses" section in
  `AGENTS.md` describing when and how to add entries (`AGENTS.md`,
  `THIRD_PARTY_NOTICES.md`)
- `Zlib_jll` added as a direct dependency (already a dependency of GZip) for
  reading unfinished gzip files (`Project.toml`, `src/GRN.jl`)

- Binary output format (`format=:bin`) — `grn()` and `makeout()` now accept
  `format=:tsv` (default) or `format=:bin`. Binary files store raw Float64
  bytes with a 16-byte header (magic `"GRNB"`, version, ncols), no text
  formatting overhead. File extensions: `.bin` (uncompressed) or `.bin.gz`
  (gzip-compressed, level 1). Binary writes use `unsafe_write` with pointer
  arithmetic — allocation-free streaming in 1 MB chunks
  (`src/io.jl`, `src/estimate.jl`)
- `read_binary(path)` — read a GRN binary output file (`.bin` or `.bin.gz`)
  into a matrix with one row per MCMC iteration. Validates the header magic
  and version; throws on corrupt files (`src/postanalysis.jl`)
- `posterior_mean`, `posterior_inclusion_probability`, `derive_indicator`,
  and `accuracy` auto-detect binary files by extension and use streaming
  binary readers — no full file load into memory (`src/postanalysis.jl`)
- `simulate(; n, m, u, seed, density, effect_sd, noise_sd, G, Q, mu, M)` —
  simulate gene expression data and genotypes from a linear gene regulatory
  network model `Y = mu + Y*G + M*Q + epsilon`. `G` is a random strictly
  upper-triangular DAG (or user-provided), `Q` is a sparse marker-gene
  effects matrix, and `epsilon` is Gaussian noise. Returns a NamedTuple
  `(Y, M, G, Q, mu)` with the simulated data and true parameters, ready to
  pass to `grn()` (`src/simulate.jl`)
- Fixed effects correction in `grn()` via the `X` keyword — an *n* × *p*
  design matrix (same for all genes) whose coefficients `Beta` (p × m) are
  estimated in the MCMC chain with flat priors. Defaults to `ones(n, 1)`
  (intercept only), so the per-gene mean is now captured as the coefficient
  of the ones column. The separate `mu`/`_sample_mean!` path is removed;
  output is always `base_fixedeffects.tsv.gz` (replaces `base_average.tsv.gz`).
  Coordinate-wise Gibbs, allocation-free
  (`src/estimate.jl`, `src/initialize.jl`, `src/sample.jl`)
- **Breaking:** `base_average.tsv.gz` is no longer written — the intercept is
  now part of `base_fixedeffects.tsv.gz`. The `mu` key in the `priors`
  NamedTuple is no longer recognized (`src/estimate.jl`, `src/initialize.jl`)
- `simulate()` gains `X` and `Beta` keywords to generate data with known
  fixed effects; returns `(Y, M, G, Q, mu, X, Beta)` (`src/simulate.jl`)

### Changed

- **Breaking:** `X` can no longer be `nothing` in `grn()` and
  `create_databucket()`. The default is now `ones(n, 1)` (intercept only)
  directly in the type signature (`X::Matrix{Float64}=ones(n, 1)`), removing
  the `Union{Matrix{Float64}, Nothing}` sentinel. The convenience method
  `grn(Y, chain; ...)` and `create_databucket()` no longer accept
  `X=nothing` (`src/estimate.jl`, `src/initialize.jl`)
- `grn(Y, chain; ...)` convenience method: `X=nothing` now defaults to
  `ones(n, 1)` (intercept only) instead of passing `nothing` to the main
  method, which expected `Matrix{Float64}` (`src/estimate.jl`)
- `AGENTS.md` rewritten in a standard structure (overview, commands, layout,
  key decisions, style, testing, definition of done, changelog, git, gotchas,
  ask-first list) and brought up to date with binary output, `simulate()`,
  and fixed effects. `/out/` added to `.gitignore` as the scratch output
  folder (`AGENTS.md`, `.gitignore`)

### Performance

- Binary output (`format=:bin`) eliminates `%.10g` text formatting and
  tab-separator overhead from the write hot path. Writes use `unsafe_write`
  with pointer arithmetic, chunked through the 1 MB IOBuffer — no per-element
  allocation (`src/io.jl`)
- `write_output` now uses a persistent 1 MB IOBuffer per output file that is
  flushed to the stream whenever it fills, instead of accumulating the entire
  row in memory and writing it in a single call. The buffer is created once
  in `makeout` and reused across all writes. Memory usage is bounded to
  ~1 MB per file regardless of row size — for 17K × 17K effects matrices
  this reduces peak memory from several GB to ~8 MB total across all output
  files (`src/io.jl`).
- Spike-and-slab effect sampler: per-gene constants `log(1-pi)`, `log(pi)`,
  and `rvar/var` are hoisted out of the per-effect loop, and the log-odds
  `logD0 - logD1` is computed directly (the `log(ssqX[q])` term cancels
  algebraically). This reduces 4 `log` calls + 1 division per effect to
  1 `log` call per effect — for 17K × 17K effects that saves ~1.7 billion
  `log` calls per MCMC iteration (`src/sample.jl`).
- `write_output` no longer flushes the IOBuffer to the stream after every
  row. Data accumulates across rows and is only written when the 1 MB
  buffer fills, batching GZip compression and reducing system call
  overhead. New `flush_outputs()` flushes pending buffer data;
  `close_outputs()` now calls `flush_outputs()` before closing streams.
  `grn()` uses `close_outputs()` in its `finally` block
  (`src/io.jl`, `src/estimate.jl`).
- Benchmark on `data/` (1305 individuals, 10 genes, 200 markers, 50K
  iterations, 981 output rows, 1 thread):

  | Configuration              | Before  | After   | Change  |
  |----------------------------|---------|---------|---------|
  | SpikeSlab (no hyperpriors) | 28.9 s  | 27.5 s  | -5%     |
  | SpikeSlab (with hyperpriors)| 30.2 s | 28.1 s  | -7%     |
  | Ridge + Laplace            | 58.0 s  | 55.9 s  | -4%     |
  | Plain output               | 29.8 s  | 26.5 s  | -11%    |

  Savings scale with problem size; the log-hoisting benefit grows
  quadratically with the number of genes (m² effects per iteration).

## [0.1.0] - 2026-10-13

Rewrite of the sampler core: the Gibbs hot path is now allocation-free, and
the prior for gene-gene and marker-gene effects is selectable per relation
type instead of being fixed to spike-and-slab.

### Added

- Prior selection with selector types in the `priors` NamedTuple
  (`src/priors.jl`):
  - `SpikeSlab(; includeprob, var, vardf)` - spike-and-slab prior (default,
    previous behavior)
  - `Ridge(; var, vardf)` - Gaussian prior on all effects with a
    sampled per-gene prior variance
  - `Laplace(; lambda2, a, b)` - Bayesian lasso prior with per-effect local
    variances (scale mixture of normals, Park and Casella 2008)
  - Selected with `priors = (Gprior = ..., Qprior = ...)`; a new prior is
    added as a selector and state subtype plus three methods
- Hyperparameters for the new priors can be set with the `priors`
  NamedTuple: `Gtau2`, `Glambda2`, `Ga`, `Gb` (and the `Q` equivalents)
- Convenience method `grn(Y, chain; ...)` for analyses without genotypes
- Marker-gene output files are no longer created when there are no marker
  loci
- Input checks for the `grn` function (`src/check.jl`), which is now
  included in the module: matching row counts, finite values, non-empty
  increasing chain, and valid prior selectors
- Test suite (`test/runtests.jl`): input checks, prior insertion, scalar
  samplers, output formats for all priors, network recovery on simulated
  data, and allocation assertions
- Benchmark scripts in `exclude/` (not part of the package): `bench.jl`
  (timing, allocations, and recovery on the data in `data/`) and
  `allocbench.jl` (per-function allocation counts)
- `posterior_inclusion_probability(path)` — compute the exact posterior
  inclusion probability (PIP) from an effects file by counting non-zero
  values across MCMC iterations (`src/postanalysis.jl`)
- `derive_indicator(path; threshold=0.5)` — derive a binary indicator
  vector from an effects file by thresholding the PIP
  (`src/postanalysis.jl`)
- `close_outputs(outs)` — close all IO streams in an output vector
  (`src/io.jl`)
- `Printf` standard library dependency for `@printf` formatting
- `AGENTS.md` — agent guidelines for the repository
- Tests for compressed output (default), plain output (`compress=false`),
  PIP derivation, and `derive_indicator` (`test/runtests.jl`)

### Changed

- `grn` runs: the main loop previously called undefined variables and
  `create_databucket` referenced undefined names, so the package did not
  run; all wiring is fixed and verified end to end
- Allocation-free Gibbs hot path:
  - Per-effect and per-hyperparameter updates are scalar loops in place
    on the data bucket; no temporary vectors are created
  - The whole per-gene sweep (mean, residual variance, gene-gene effects,
    marker-gene effects, and hyperparameters) is fused into a single
    threaded region per iteration, since each stage depends only on the
    residual of its own gene
  - Chi-squared, beta, and gamma draws use a scalar Marsaglia and Tsang
    (2000) gamma sampler; the Laplace local variances use the Michael,
    Schucany, and Haas (1976) inverse Gaussian sampler
  - The remaining allocation per iteration is the task creation of
    `@threads` (about 0.9 kB per task); a persistent task pool would be
    needed to remove it
  - Measured on `data/` (1305 individuals, 10 genes, 200 markers,
    12 threads): 0 bytes for the Gibbs updates, about 8 kB per iteration
    total, 0.28-0.41 ms per iteration
- The sweep order changed from staged (all means, then all residual
  variances, then all gene-gene effects, then all marker-gene effects) to
  per-gene fused. Both are valid Gibbs orderings, but results are not
  identical to previous runs even with the same seed
- `insert_priors!` maps keys to fields by name (base fields by field name,
  prior fields prefixed with the group letter) and recalculates the prior
  variance scales unless a `varscale` key is given explicitly
- Output files are created from a flat list of output objects: one
  `<group>_<parameter>.tsv` file per object, with names that depend on the
  selected prior (for example `G_tau2.tsv` and `G_lambda2.tsv` only for the
  Laplace prior)
- Residual variance defaults: the prior residual variance scale is now
  initialized consistently as `Rvar * (Rvardf - 2)`
- **Breaking:** `compress=true` is now the default in `grn()` and
  `makeout()`. Output files are gzip-compressed (`.tsv.gz`, level 1) by
  default. Set `compress=false` for plain `.tsv` files (`src/estimate.jl`,
  `src/io.jl`)
- File writes use an `IOBuffer` with `sizehint` pre-allocation and
  `@printf("%.10g", ...)` for Float64 formatting, replacing per-element
  `print` calls to the IO stream. `_write_row!` methods are specialized
  for `Float64` and nested `Float64` vectors (`src/io.jl`)
- `GZip.open` uses positional compression level argument (`"w", 1`)
  matching GZip.jl v0.8 API (`src/io.jl`)
- `accuracy`, `posterior_mean`, `posterior_inclusion_probability`, and
  `derive_indicator` are now exported from the module
- Test helper `readoutput` handles both `.tsv` and `.tsv.gz` files
  (`test/runtests.jl`)

### Removed

- `Distributions` dependency: all draws are generated with scalar samplers
  internal to the package
- Internal helpers that are no longer used: `_sample_pi!`,
  `_sample_effectvars!`, `dot!`, `indicator2indice`, and the temporary
  per-thread buffers (`GFloat*`, `QFloat*`, `IFloat*`) from the data bucket
- **Breaking:** Indicator files (`G_indicator.tsv.gz`, `Q_indicator.tsv.gz`)
  are no longer written. The posterior inclusion probability is derived
  exactly from the effects files via `posterior_inclusion_probability()`
  (`src/priors.jl`, `src/estimate.jl`)
- Old `_write_row` methods (non-buffer, direct-to-stream) — replaced by
  `_write_row!` (IOBuffer-based) in `src/io.jl`
- Int8 fast path in `_write_row!` — dead code after indicator removal

### Performance

Comparison on the simulated data in `data/` (10000 iterations with 5000
burn-in, posterior mean over 51 kept samples, recovery measured as the
correlation with the true effects in `B.csv` and `F.csv`):

| Prior      | Time  | cor(G) | cor(Q) |
|------------|-------|--------|--------|
| SpikeSlab  | 5.9 s | 0.770  | 0.468  |
| Laplace    | 12.4 s | 0.805  | 0.518  |
| Ridge      | 10.2 s | 0.739  | 0.444  |

- Faster file writes: `@printf` Float64 formatting + `IOBuffer` batching
  with `sizehint` pre-allocation replaces per-element `print` calls to the
  IO stream (`src/io.jl`). For a 17K x 17K effects matrix (~289M values per
  row), the row is buffered in memory and written in a single `write` call.
