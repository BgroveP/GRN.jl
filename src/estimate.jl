"""

    grn(Y, M, chain; X=ones(size(Y,1), 1), priors=(;), outpath="grn", overwrite=false)

Estimate the gene regulatory network using a matrix with gene expression and
a matrix with genotypes.


# Dimensions
 - *n*:           The number of individuals
 - *m*:           The number of genes
 - *u*:           The number of marker loci


# Mandatory arguments
 - `Y`:             An *n* × *m* matrix with gene expressions                                                           (Float64)
 - `M`:             An *n* × *u* matrix with genotypes                                                                  (Float64)

# Optional arguments
 - `chain`:         The MCMC iterations where parameters are written to the output files                                (AbstractRange{Int})
                    For example, giving `10000:50:20000` as input will regard iterations
                       1:9999 as burn-in, and write parameters for iterations 10000,
                       10050, 10100, ..., 20000
 - `X`:             An *n* × *p* matrix with fixed effect covariates, the same
                    for all genes. Defaults to `ones(n, 1)` (intercept only).
                    The coefficients `Beta` (p × m) are estimated in the MCMC
                    chain with flat priors. The per-gene intercept is the
                    coefficient of the ones column.                                                     (Matrix{Float64})
 - `priors`:        The priors given for the analysis                                                                   (NamedTuple)
 - `outpath`:       The relative path to the folder for output files                                                    (String)
 - `overwrite`:     Denotes whether an existing outpath should be overwritten                                           (Bool)
 - `compress`:      Gzip-compress the output files (default, level 1, fast). Set to
                    `false` for uncompressed files.                                                        (Bool)
 - `format`:        Output format: `:tsv` for tab-separated text (default) or `:bin` for
                    raw binary Float64 with a 16-byte header. Binary format is much
                    faster to write and read — recommended for large analyses.
                    File extensions: `.tsv`, `.tsv.gz`, `.bin`, `.bin.gz`.                                    (Symbol)

# Priors
 The prior for gene-gene (G) and marker-gene (Q) effects is selected with
 the `Gprior` and `Qprior` keys, which take a prior selector:

  - `SpikeSlab(; includeprob, var, vardf)`: Spike-and-slab prior with sampled
    inclusion probabilities (default)
  - `Ridge(; var, vardf)`: Gaussian prior on all effects with a sampled
    per-gene prior variance
  - `Laplace(; lambda2, a, b)`: Bayesian lasso prior with per-effect local
    variances and a sampled per-gene penalty `lambda2`

 For example, `priors = (Gprior = Laplace(), Qprior = Ridge())`. All other
 keys set hyperparameters of the model, and are inserted into the state
 before the analysis. The key for a base field is the field name, and the
 key for a prior state field is the group letter followed by the field name.

 Priors related to the residual
 - `Rvar`:          An *m* dimensional vector of residual variances                                                     (Float64)
 - `Rvardf`:        An *m* dimensional vector of degrees of freedoms for the prior residual variance                    (Float64)
 - `Rvarscale`:     An *m* dimensional vector of prior degrees of freedom times prior residual variances                 (Float64)

 Priors related to gene-gene effects (replace `G` with `Q` for marker-gene effects)
 - `Geffects`:      An *m* × *m* dimensional matrix of effects                                                          (Float64)
 - `Gpi`:          An *m* dimensional vector of the probability that an effect is not zero                             (Float64)
 - `Gvar`:          An *m* dimensional vector of variances for effects                                                  (Float64)
 - `Gvardf`:        An *m* dimensional vector of degrees of freedom for the prior variance                              (Float64)
 - `Gvarscale`:     An *m* dimensional vector of prior degrees of freedom times prior effect variances                  (Float64)
 - `Gtau2`:         An *m* dimensional vector of **m** dimensional vectors of local variances (Laplace)                 (Float64)
 - `Glambda2`:      An *m* dimensional vector of squared penalties (Laplace)                                            (Float64)
 - `Ga`:            An *m* dimensional vector of shapes for the prior of lambda2 (Laplace)                               (Float64)
 - `Gb`:            An *m* dimensional vector of rates for the prior of lambda2 (Laplace)                                (Float64)
 - `Gcandidates`:   An *m* dimensional vector of vectors with indices of candidate genes                                (Int)

# Continuing an analysis
 `grn(path, Y, M, chain; ...)`, with the output folder of an earlier analysis
 as first argument, continues that analysis from its last stored iteration.
 See the docstring of that method.

# Returns
 - `Nothing`: The relevant output is written to files

"""
function grn(Y::Matrix{Float64},
    M::Matrix{Float64}=zeros(0, 0),
    chain::AbstractRange{Int}=100:10:200;
    X::Matrix{Float64}=ones(size(Y, 1), 1),
    priors::NamedTuple=(;),
    outpath::AbstractString="grn",
    overwrite::Bool=false,
    compress::Bool=true,
    format::Symbol=:tsv)

    # Check input
    _checkinput(Y, M, chain, X)

    # Initialize the data bucket with priors inserted
    data, amounts, ssqs = _setup(Y, M, X, priors)

    # Create folder and IO streams
    ## This function creates a named collection of IOStreams for later use
    ofiles = makeout(outpath, _outputs(data, amounts.loci); force=overwrite, compress=compress, format=format)

    # Run the MCMC chain
    _run_chain!(data, ssqs, amounts, chain, ofiles)

    # Return
    return nothing
end

"""

    _setup(Y, M, X, priors)

Print the dimensions, select the priors, create the data bucket, insert the
priors, and compute the squared sums of the predictors. Shared by new and
continued analyses. Returns `(data, amounts, ssqs)`.

"""
function _setup(Y::Matrix{Float64}, M::Matrix{Float64}, X::Matrix{Float64},
    priors::NamedTuple)::Tuple{NamedTuple, NamedTuple, NamedTuple}

    # Get amounts of individuals, genes, marker loci, and fixed effects
    amounts = (individuals=size(Y, 1), genes=size(Y, 2), loci=size(M, 2),
        fixed=size(X, 2))

    # Print
    _print(amounts)

    # Select the priors for gene-gene and marker-gene effects
    Gsel = get(priors, :Gprior, SpikeSlab())
    Qsel = get(priors, :Qprior, SpikeSlab())
    Gsel isa AbstractEffectPrior || throw(ArgumentError(
        "Gprior must be a prior selector (SpikeSlab(), Ridge(), or Laplace()), got $(typeof(Gsel))"))
    Qsel isa AbstractEffectPrior || throw(ArgumentError(
        "Qprior must be a prior selector (SpikeSlab(), Ridge(), or Laplace()), got $(typeof(Qsel))"))

    # Initialize parameters
    ## This creates a large named tuple with input, output, and temporary arrays
    data = create_databucket(Y, M, amounts, Gsel, Qsel, X)

    # Insert priors
    ## This function inserts priors into the output part of the data bucket
    insert_priors!(data, priors)

    # Calculate squared sums for each column of Y and M
    ssqs = (
        M=[dot(m, m) for m in data.in.M],
        Y=[dot(y, y) for y in data.in.Y]
    )

    return data, amounts, ssqs
end

"""

    _run_chain!(data, ssqs, amounts, chain, ofiles)

Run `last(chain)` Gibbs iterations and write the outputs for the iterations
in `chain`. The output files are closed afterwards, also on error.

"""
function _run_chain!(data::NamedTuple, ssqs::NamedTuple, amounts::NamedTuple,
    chain::AbstractRange{Int}, ofiles::Vector{NamedTuple})::Nothing

    try
        # Start analysis
        for i in 1:last(chain)

            # Sample parameters
            ## The full Gibbs sweep is performed per gene, in parallel
            ## over contiguous chunks of genes
            _sample_iteration!(data, ssqs, amounts.individuals, amounts.loci)

            # Write to file if iteration is an output iteration
            if i in chain
                @threads for k in eachindex(ofiles)
                    write_output(ofiles[k])
                end
            end
        end
    finally
        # Flush buffers and close the output files
        close_outputs(ofiles)
    end

    return nothing
end

# Convenience method for analyses without genotypes
function grn(Y::Matrix{Float64}, chain::AbstractRange{Int};
    X::Matrix{Float64}=ones(size(Y, 1), 1),
    priors::NamedTuple=(;),
    outpath::AbstractString="grn",
    overwrite::Bool=false,
    compress::Bool=true,
    format::Symbol=:tsv)
    return grn(Y, zeros(0, 0), chain; X=X, priors=priors, outpath=outpath,
        overwrite=overwrite, compress=compress, format=format)
end

"""

    grn(path, Y, M, chain; X=ones(size(Y,1), 1), priors=(;), outpath=path, overwrite=false)

Continue an earlier analysis from the results stored in the folder `path`.
This is the method of `grn` used when the first argument is a folder
instead of the expression matrix. The sampler state is restored from the last iteration stored in the
output files, and `last(chain)` further iterations are run. By default the
new iterations are appended to the existing files, so the files hold one
long chain that the post-analysis functions read as usual.

The output files hold the sampled parameters, but not the data or the fixed
prior settings. Give the same `Y`, `M`, `X`, and `priors` as in the original
call to `grn`, so that the continued chain samples the same model. The prior
types (`Gprior`, `Qprior`) and the output format and compression are taken
from the files; a `Gprior` or `Qprior` in `priors` must match them.

# Arguments
 - `path`:          The folder with the output files of the earlier analysis                                            (String)
 - `Y`, `M`, `X`:   The data of the earlier analysis (see `grn`)                                                        (Matrix{Float64})
 - `chain`:         The iterations to write, counted from the restart. For example, `10:10:1000` runs
                    1000 more iterations and writes every 10th. A burn-in is usually not needed, since the
                    chain continues from a sample of the earlier run.                                                   (AbstractRange{Int})
 - `priors`:        The priors of the earlier analysis (see `grn`)                                                      (NamedTuple)
 - `outpath`:       Where to write the new iterations. With the default (`path`), they are appended to the
                    existing files. Any other folder gets new files with only the new iterations.                       (String)
 - `overwrite`:     Denotes whether an existing `outpath` other than `path` should be overwritten                      (Bool)

# Interrupted analyses
If the earlier analysis was stopped before it finished (e.g. by a time
limit), the files can end with an incomplete row or an unfinished gzip
stream, and the files can hold different numbers of iterations, since each
file is buffered separately. The state is then restored from the last
iteration that is complete in all files. Appending to such files is not
possible, so give a new `outpath`.

# Notes
 - Text output stores values with 10 significant digits, so the restored
   state is rounded to that precision. Binary output (`format=:bin`) is
   restored exactly.
 - The random number generator state is not stored. The continued chain is
   a valid continuation of the Markov chain, but not the same draws as one
   uninterrupted run.
 - Restoring reads every output file once (twice for files with more
   iterations than the others after an interruption).

# Returns
 - `Nothing`: The relevant output is written to files

# Example
```julia
grn(Y, M, 10000:50:100000; outpath="grn_example")
grn("grn_example", Y, M, 50:50:50000)    # 1000 more stored iterations
```

"""
function grn(path::AbstractString, Y::Matrix{Float64},
    M::Matrix{Float64}=zeros(0, 0),
    chain::AbstractRange{Int}=10:10:100;
    X::Matrix{Float64}=ones(size(Y, 1), 1),
    priors::NamedTuple=(;),
    outpath::AbstractString=path,
    overwrite::Bool=false)

    # Check input
    isdir(path) || throw(ArgumentError("No output directory at $path"))
    _checkinput(Y, M, chain, X)
    append = ispath(outpath) && samefile(outpath, path)
    overwrite && _inside(path, outpath) && throw(ArgumentError(
        "overwrite=true would delete the output files in $path that are being continued"))

    # Restore the state of the last stored iteration
    st = _continued_state(path, Y, M, X, priors)

    # Appending needs files that were closed properly
    if append && !st.complete
        error("The output files in $path were not closed properly (the analysis was " *
              "probably interrupted), so new iterations cannot be appended to them. " *
              "Give a new folder with `outpath` to continue from iteration $(st.iteration).")
    end

    # Open the output files and run the chain
    ofiles = makeout(outpath, _outputs(st.data, st.amounts.loci); force=overwrite,
        compress=st.compress, format=st.format, append=append)
    _run_chain!(st.data, st.ssqs, st.amounts, chain, ofiles)

    return nothing
end

# Convenience method for analyses without genotypes
function grn(path::AbstractString, Y::Matrix{Float64}, chain::AbstractRange{Int};
    X::Matrix{Float64}=ones(size(Y, 1), 1),
    priors::NamedTuple=(;),
    outpath::AbstractString=path,
    overwrite::Bool=false)
    return grn(path, Y, zeros(0, 0), chain; X=X, priors=priors,
        outpath=outpath, overwrite=overwrite)
end

"""

    _inside(inner, outer)

Return `true` if the existing path `inner` is `outer` or lies inside it.

"""
function _inside(inner::AbstractString, outer::AbstractString)::Bool
    ispath(outer) || return false
    p = realpath(inner)
    while true
        samefile(p, outer) && return true
        parent = dirname(p)
        parent == p && return false
        p = parent
    end
end

"""

    _continued_state(path, Y, M, X, priors)

Set up the data bucket as `grn` does, then restore the sampled parameters
from the last iteration stored in all output files in `path`, rebuild the
state that is not written to file (`_restore!`), and recompute the
residuals. Returns a NamedTuple with `data`, `amounts`, `ssqs`, the
`format` and `compress` of the files, the restored `iteration` (row number
in the files), and `complete` (all files ended cleanly with the same number
of iterations).

"""
function _continued_state(path::AbstractString, Y::Matrix{Float64},
    M::Matrix{Float64}, X::Matrix{Float64}, priors::NamedTuple)::NamedTuple

    # Format and compression of the earlier analysis
    ext = _output_extension(path)
    format = startswith(ext, ".bin") ? :bin : :tsv
    compress = endswith(ext, ".gz")

    # Prior types of the earlier analysis
    priors = _continued_priors(path, ext, priors)

    # Set up the data bucket exactly as grn does
    data, amounts, ssqs = _setup(Y, M, X, priors)
    outputs = _outputs(data, amounts.loci)
    _check_output_files(path, ext, keys(outputs))

    # Restore the sampled parameters and the state derived from them
    iteration, complete = _restore_outputs!(outputs, path, ext)
    _restore!(data.out.G)
    amounts.loci > 0 && _restore!(data.out.Q)
    _residuals!(data, Y, amounts.loci)
    println("Continuing from stored iteration $iteration in $path")

    return (data=data, amounts=amounts, ssqs=ssqs, format=format,
        compress=compress, iteration=iteration, complete=complete)
end

"""

    _output_extension(path)

Return the file extension (`.tsv.gz`, `.tsv`, `.bin.gz`, or `.bin`) of the
output files in `path`, from the `base_residual_variance` file that every
analysis writes.

"""
function _output_extension(path::AbstractString)::String
    found = [e for e in (".tsv.gz", ".tsv", ".bin.gz", ".bin")
             if isfile(joinpath(path, "base_residual_variance$e"))]
    isempty(found) && throw(ArgumentError(
        "No GRN output files in $path (no base_residual_variance file)"))
    length(found) == 1 || throw(ArgumentError(
        "Output files in several formats in $path: $(join(found, ", "))"))
    return found[1]
end

"""

    _continued_priors(path, ext, priors)

Add the prior selectors of the earlier analysis in `path` to `priors`, or
check that the selectors given in `priors` are of the same type.

"""
function _continued_priors(path::AbstractString, ext::AbstractString,
    priors::NamedTuple)::NamedTuple
    for (g, key) in ((:G, :Gprior), (:Q, :Qprior))
        T = _detect_prior(path, g, ext)
        T === nothing && continue
        if haskey(priors, key)
            priors[key] isa T || throw(ArgumentError(
                "$key is $(typeof(priors[key])), but the analysis in $path used $T"))
        else
            priors = merge(priors, NamedTuple{(key,)}((T(),)))
        end
    end
    return priors
end

"""

    _check_output_files(path, ext, names)

Check that `path` holds exactly the output files `names` (with extension
`ext`) among its `base_`, `G_`, and `Q_` files, so the data and priors match
the earlier analysis.

"""
function _check_output_files(path::AbstractString, ext::AbstractString, names)::Nothing
    expected = ["$k$ext" for k in names]
    for f in expected
        isfile(joinpath(path, f)) || throw(ArgumentError(
            "Missing output file $f in $path. Y, M, X, and priors must match the analysis being continued"))
    end
    extra = [f for f in readdir(path) if endswith(f, ext) && !(f in expected) &&
             any(p -> startswith(f, p), ("base_", "G_", "Q_"))]
    isempty(extra) || throw(ArgumentError(
        "Output files $(join(extra, ", ")) in $path are not written by an analysis with the " *
        "given data and priors (e.g. Q_ files need the genotypes M)"))
    return nothing
end

"""

    _restore_outputs!(outputs, path, ext)

Overwrite each output object in `outputs` with the last iteration that is
stored in all output files in `path`. Files are read in parallel, and files
holding more iterations are read again up to that iteration. Returns the
iteration (row number) and whether all files ended cleanly with the same
number of iterations.

"""
function _restore_outputs!(outputs::NamedTuple, path::AbstractString,
    ext::AbstractString)::Tuple{Int, Bool}
    names = collect(keys(outputs))
    files = [joinpath(path, "$k$ext") for k in names]
    nrows = zeros(Int, length(names))
    clean = [Ref(false) for _ in names]

    try
        # Last row of every file
        @threads for j in eachindex(names)
            nrows[j] = _read_row!(outputs[names[j]], files[j], 0; complete=clean[j])
        end
        minimum(nrows) > 0 || error("No iteration is stored in all output files in $path")

        # Files that are ahead of the others are read up to iteration n
        n = minimum(nrows)
        @threads for j in findall(>(n), nrows)
            _read_row!(outputs[names[j]], files[j], n)
        end
    catch e
        _rethrow_task_error(e)
    end

    n = minimum(nrows)
    return n, all(r -> r[], clean) && all(==(n), nrows)
end

"""

    _rethrow_task_error(e)

Rethrow the original exception from inside a `@threads` loop, which arrives
wrapped in a `TaskFailedException` or `CompositeException`.

"""
function _rethrow_task_error(e)
    while true
        if e isa CompositeException && !isempty(e.exceptions)
            e = first(e.exceptions)
        elseif e isa TaskFailedException
            e = e.task.exception
        else
            throw(e)
        end
    end
end

"""

    _read_row!(dest, file, target; complete=nothing)

Copy row `target` of an output file (the last complete row for
`target == 0`) into the output object `dest`, and return the number of
complete rows in the file. The file is streamed, so memory use is one row.

"""
function _read_row!(dest::AbstractVector, file::AbstractString, target::Int;
    complete::Union{Nothing, Base.RefValue{Bool}}=nothing)::Int
    ncols = _ncols(dest)
    i = Ref(0)
    return _scan_rows(file, '\t', nothing; complete=complete) do row
        length(row) == ncols || throw(DimensionMismatch(
            "$file has $(length(row)) values per iteration, but $ncols are expected " *
            "for the given data and priors"))
        i[] += 1
        (target == 0 || i[] == target) && _copyrow!(dest, row)
        return nothing
    end
end

"""

    _copyrow!(dest, row)

Copy a flat row of values into an output object: a vector, or a vector of
vectors filled in order (the layout `write_output` uses).

"""
function _copyrow!(dest::Vector{Float64}, row::AbstractVector{Float64})::Nothing
    copyto!(dest, row)
    return nothing
end

function _copyrow!(dest::Vector{Vector{Float64}}, row::AbstractVector{Float64})::Nothing
    o = 0
    for v in dest
        copyto!(v, 1, row, o + 1, length(v))
        o += length(v)
    end
    return nothing
end

"""

    _residuals!(d, Y, loci)

Recompute the residual of every gene from the data and the current
parameters: `y - X * beta - Yc * G - Mc * Q`, with `Yc` and `Mc` the
centered predictors used by the sampler. Run in parallel over genes; each
task writes only the residuals of its own genes.

"""
function _residuals!(d::NamedTuple, Y::Matrix{Float64}, loci::Int)::Nothing
    chunks = d.tmpshared.split_chunks
    @threads for s in eachindex(chunks)
        for c in chunks[s]
            r = d.in.R[c]
            r .= view(Y, :, c)
            _subtract_effects!(r, d.out.base.fixedeffects[c], d.in.X)
            _subtract_effects!(r, d.out.G.effects[c], d.in.Y)
            loci > 0 && _subtract_effects!(r, d.out.Q.effects[c], d.in.M)
        end
    end
    return nothing
end

"""

    _subtract_effects!(r, w, X)

Subtract `X[q] * w[q]` from `r` for every non-zero effect `w[q]`.

"""
function _subtract_effects!(r::Vector{Float64}, w::Vector{Float64},
    X::Vector{Vector{Float64}})::Nothing
    @inbounds for q in eachindex(w)
        wq = w[q]
        wq != 0.0 && (r .-= wq .* X[q])
    end
    return nothing
end

"""

    _outputs(d, loci)

Collect the output objects of the data bucket into a flat NamedTuple with
one entry per output file. The keys become the file names. The marker-gene
group is skipped when there are no marker loci. The fixed effects
(`base_fixedeffects`) include the intercept when `X` is a column of ones.
Indicator matrices are not included — the posterior inclusion probability
can be derived exactly from the effects files with
`GRN.posterior_inclusion_probability`.

"""
function _outputs(d::NamedTuple, loci::Int)::NamedTuple

    # Base outputs
    o = (
        base_fixedeffects=d.out.base.fixedeffects,
        base_residual_variance=d.out.base.residual_variance
    )

    # Prior state outputs
    groups = loci > 0 ? (:G, :Q) : (:G,)
    for g in groups
        st = d.out[g]
        so = _outputs(st)
        for k in keys(so)
            o = merge(o, NamedTuple{(Symbol(g, "_", k),)}((so[k],)))
        end
    end

    return o
end
