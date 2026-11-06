# AGENTS.md — GRN.jl

Instructions for AI coding agents (Claude Code, Mistral Vibe, and others)
working in this repository. Humans are welcome to read it too.

## Project overview

GRN.jl is a Julia package for Bayesian gene regulatory network inference by
Gibbs sampling. Each gene's expression is regressed on the expression of the
other genes (gene-gene effects, `G`), optionally on marker genotypes
(marker-gene effects, `Q`), and on fixed-effect covariates `X`. The priors on
`G` and `Q` are selectable (`SpikeSlab`, `Ridge`, `Laplace`). Sampled
parameters are streamed to output files, one row per retained iteration.

**Scale target:** 17K genes × 17K markers. One effects row can be ~5 GB of
uncompressed text per iteration. Memory use and write speed matter as much as
sampler speed. Profile before optimizing.

## Commands

Run from the package root. The developer works on Windows (Git Bash or
PowerShell, Julia via juliaup); CI runs Ubuntu on Julia 1.10, 1.12, and pre.

```sh
julia --project -e 'using Pkg; Pkg.instantiate()'          # set up
julia --project -e 'using Pkg; Pkg.test()'                 # full test suite
julia --project -t 4 -e 'using Pkg; Pkg.test()'            # exercise the threaded path
julia --project examples/examples.jl                       # end-to-end examples (slow, writes grn_example*/)
julia --project=docs docs/make.jl                          # build docs
```

- Code must run on **Julia 1.10** (the `[compat]` lower bound). Do not use
  syntax or Base functions introduced later.
- `Manifest.toml` is gitignored. Do not commit it.

## Repository layout

| Path | Contents |
|---|---|
| `src/GRN.jl` | Module definition, imports, exports, `_BUFSIZE` (1 MB write buffer) |
| `src/priors.jl` | Prior selectors (`SpikeSlab`, `Ridge`, `Laplace`), state types, `_initstate`, `_outputs`, `_varscale!` |
| `src/initialize.jl` | `create_databucket`, `insert_priors!` |
| `src/sample.jl` | Gibbs samplers and internal scalar random draws (allocation-free hot path) |
| `src/estimate.jl` | `grn()` entry point, `_outputs()` collector |
| `src/io.jl` | Output files: `makeout`, `write_output`, `flush_outputs`, `close_outputs`; TSV and binary writers |
| `src/postanalysis.jl` | `accuracy`, `posterior_mean`, `posterior_inclusion_probability`, `derive_indicator`, `read_binary`, `preliminary_posterior` (files still being written) |
| `src/simulate.jl` | `simulate()` — data from a known network, for tests and validation |
| `src/check.jl`, `src/utility.jl` | Input validation, helpers |
| `test/runtests.jl` | Single test file, grouped in `@testset`s |
| `examples/examples.jl` | End-to-end runs on `data/` with accuracy report |
| `data/` | Example input (`Y.csv`, `M.csv`) and true effects (`B.csv` for G, `F.csv` for Q) |
| `docs/` | Documenter.jl site |

## Key decisions — do not regress

- **No indicator files are written.** The posterior inclusion probability is
  derived exactly from the effects files by `posterior_inclusion_probability()`
  (share of non-zero values per column), so indicator files would only add
  output volume.
- **Output is gzip-compressed by default** (`compress=true`, level 1) in both
  `grn()` and `makeout()`. `compress=false` gives plain files.
- **Two output formats:** `format=:tsv` (default) and `format=:bin`. Binary
  files have a 16-byte header (magic `"GRNB"`, version `UInt32`, ncols
  `Int64`) followed by raw native-endian `Float64` rows. Extensions:
  `.tsv`, `.tsv.gz`, `.bin`, `.bin.gz`. Any change to the binary layout must
  bump `_BIN_VERSION` and keep `read_binary` able to reject unknown versions.
- **Writes go through a persistent 1 MB `IOBuffer` per file**, flushed only
  when full, so memory stays bounded regardless of row size. Text uses
  `@printf` with `%.10g` for `Float64`. Never `print` element by element
  directly to a stream.
- **Post-analysis functions stream files.** They must not load a whole output
  file into memory.
- **The Gibbs hot path is allocation-free.** Scalar loops, in-place updates,
  no temporary vectors. This is enforced by the `Allocations` testset
  (`@allocated == 0` for the per-gene updates). Cold paths (setup, generic
  fallbacks for non-`Float64` types) may allocate.
- **No `Distributions` dependency.** All random draws use the internal scalar
  samplers in `src/sample.jl` (`_rand_gamma`, `_rand_beta`, `_rand_invgauss`).
- **Fixed effects:** `X` defaults to `ones(n, 1)`. The intercept is the
  coefficient of the ones column; there is no separate mean parameter.
  `X` is never `nothing`.

## Adding a new prior

1. In `src/priors.jl`: define a selector `<: AbstractEffectPrior` (with a
   keyword constructor and docstring) and a state `<: AbstractEffectPriorState`.
2. Add methods for `_initstate`, `_outputs`, and `_varscale!` (`src/priors.jl`)
   and `_sample_effects_gene!`, `_sample_hyperparameters_gene!` (`src/sample.jl`).
3. Export the selector in `src/GRN.jl`, list it in the `grn` docstring, README,
   and `docs/src/index.md`.
4. Add tests: output files written, allocation-free updates, effect recovery
   on `simulate()` data.

## Code style

- 4-space indentation, no tabs.
- Type annotations on function signatures and return types (`::Float64`,
  `::Nothing`, …).
- `@inbounds` in hot loops.
- `"""` docstrings for all public **and** internal functions.
- Immutable structs and NamedTuples for state; no mutable structs. Mutate
  the contents of arrays in place instead.
- Follow the existing comment style: short `# Section` comments above blocks.

## Testing

- All tests live in `test/runtests.jl`. Add a new `@testset` or extend the
  matching one; do not create extra test files without asking.
- Current testsets: input checks, prior selection, scalar samplers, simulate,
  fixed effects, prior insertion, analysis and output, compressed output,
  plain output, PIP and `derive_indicator`, binary output, preliminary
  results, allocations.
- Every new output file, output format, or post-analysis function needs a
  test.
- Changes to the samplers must keep the effect-recovery tests (correlation
  with `simulate()` truth) passing. Do not loosen a recovery threshold or an
  allocation limit to make a test pass — report it instead.
- Tests write into a temporary directory; never into the repo.

## Definition of done

Before reporting a task as finished:

1. `Pkg.test()` passes. If threading code changed, it also passes with `-t 4`.
2. The `Allocations` testset still passes unchanged.
3. `CHANGELOG.md` is updated (see below).
4. If the public API or output files changed: the `grn` docstring,
   `README.md`, and `docs/src/index.md` are updated in the same change, and
   `examples/examples.jl` still runs.
5. Performance claims include before/after numbers with the setup used
   (data size, iterations, threads, Julia version).

If a step was skipped or fails, say so explicitly.

## Changelog

- **Every change must be recorded in `CHANGELOG.md`**, in the same commit as
  the code change — never in a separate commit.
- Use the existing `[Unreleased]` section at the top. Do not create a version
  header unless releasing.
- Before adding an entry, read `[Unreleased]`. If an existing entry describes
  the same feature, edit or replace it instead of adding a contradictory one.
- Subsections:
  - `### Added` — new features, functions, files, priors, options
  - `### Changed` — behavior changes, refactors, signature changes, defaults
  - `### Removed` — deleted code, dependencies, output files
  - `### Performance` — speed or memory improvements (with before/after numbers)
- One bullet per change, naming the function and file, e.g.
  `` `grn()` now accepts `compress=true` (`src/estimate.jl`) ``.
- Breaking changes start with `**Breaking:**` and go under the subsection
  that matches what happened (usually Changed or Removed), not Added.

## Dependencies

- Current deps: `Printf`, `LinearAlgebra`, `Random`, `Statistics`,
  `DelimitedFiles`, `GZip`, `Zlib_jll` (raw `gzread` for unfinished gzip
  files, see `preliminary_posterior`).
- Ask before adding a dependency. When adding one, update `Project.toml`
  `[deps]` and `[compat]`, then run `Pkg.resolve()`.

## Git

- Work on a branch, not `main`, unless told otherwise. Commit only when asked.
- One logical change per commit (e.g. a new output format and a sampler
  speed-up are two commits), each with its own CHANGELOG entry.
- Never commit output folders (`grn_example*/`, `out/`) or `Manifest.toml`.

## Gotchas

- `BLAS.set_num_threads(1)` is called at module load — the package
  parallelizes over genes with `Threads.@threads`, not through BLAS.
- `_sample_iteration!` runs the full per-gene sweep in parallel over
  contiguous chunks of genes. Code inside it may only write to data owned by
  gene `c` (its residual, its effects, its hyperparameters). Shared scratch
  buffers would race.
- Results are not bit-reproducible across different thread counts.
- `makeout(...; force=true)` deletes the whole output directory. Never point
  `outpath` at a directory that holds anything else.
- The repo lives in a OneDrive-synced folder. Write scratch and analysis
  output to `out/` or a directory outside the repo, and avoid writing large
  files into the working tree.
- `data/B.csv` is the true gene-gene matrix for `data/Y.csv` and
  `data/F.csv` the true marker-gene matrix. `accuracy(truth, effectsfile)`
  reshapes the posterior mean to the shape of the truth file and returns
  their correlation; it does not transpose.

## Ask first

- Changing a default (`compress`, `format`, prior hyperparameters, `X`).
- Changing an output file name, the binary header, or the TSV row layout.
- Adding a dependency, a new public export, or a new top-level directory.
- Anything that makes a hot-path function allocate.
