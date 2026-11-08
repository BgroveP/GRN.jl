
# Load packages
using DelimitedFiles
using GRN
using GZip

# Load the input
pkgroot = dirname(dirname(pathof(GRN)))
Y = readdlm("$(pkgroot)/data/Y.csv", ',');
M = readdlm("$(pkgroot)/data/M.csv", ',');

# Iterations 
iterations = 1000:10:3000
# Spike-and-slab priors, no other priors (compressed by default)
@time grn(Y, M, 10000:50:100000, outpath = "grn_example1", overwrite = true)

# Spike-and-slab priors, with hyperparameter priors
@time grn(Y, M, 10000:50:100000, outpath = "grn_example2", overwrite = true, priors = (Gpi = 0.01, Qpi = 0.01, Rvar = 0.01))

# Bayesian ridge prior for the gene-gene effects, lasso prior for the marker-gene effects
@time grn(Y, M, 10000:50:100000, outpath = "grn_example3", overwrite = true, priors = (Gprior = Ridge(), Qprior = Laplace()))

# Pure Ridge regression
@time grn(Y, M, 10000:50:100000, outpath = "grn_example4", overwrite = true, priors = (Gprior = Ridge(), Qprior = Ridge()))

# Writing to uncompressed plain text files
@time grn(Y, M, 10000:50:100000, outpath = "grn_example5", overwrite = true, compress = false)

# Writing to binary files
@time grn(Y, M, 10000:50:100000, outpath = "grn_example6", overwrite = true, compress = false, format=:bin)

# Writing to binary files
@time grn(Y, M, 10000:50:100000, outpath = "grn_example7", overwrite = true, format=:bin)

# Continuing grn_example1 from its last stored iteration, with one more
# iteration written to a new folder (grn_example1 itself is left unchanged)
@time grn("grn_example1", Y, M, 1:1:1, outpath = "grn_example8", overwrite = true)

# Accuracies
@info "
    Accuracies:
        Gene effects:
            Spike and slab (no priors)                      $(round(GRN.accuracy("$(pkgroot)/data/B.csv", "grn_example1/G_effects.tsv.gz"), digits = 2))
            Spike and slab (with priors)                    $(round(GRN.accuracy("$(pkgroot)/data/B.csv", "grn_example2/G_effects.tsv.gz"), digits = 2))
            Lasso+                                          $(round(GRN.accuracy("$(pkgroot)/data/B.csv", "grn_example3/G_effects.tsv.gz"), digits = 2))
            Ridge                                           $(round(GRN.accuracy("$(pkgroot)/data/B.csv", "grn_example4/G_effects.tsv.gz"), digits = 2))
            Spike and slab (no priors, plain)               $(round(GRN.accuracy("$(pkgroot)/data/B.csv", "grn_example5/G_effects.tsv"), digits = 2))
            Spike and slab (no priors, binary)              $(round(GRN.accuracy("$(pkgroot)/data/B.csv", "grn_example6/G_effects.bin"), digits = 2))
            Spike and slab (no priors, compr. binary)       $(round(GRN.accuracy("$(pkgroot)/data/B.csv", "grn_example7/G_effects.bin.gz"), digits = 2))
            Spike and slab (continued from example 1)       $(round(GRN.accuracy("$(pkgroot)/data/B.csv", "grn_example8/G_effects.tsv.gz"), digits = 2))

        Marker effects:
            Spike and slab (no priors)                      $(round(GRN.accuracy("$(pkgroot)/data/F.csv", "grn_example1/Q_effects.tsv.gz"), digits = 2))
            Spike and slab (with priors)                    $(round(GRN.accuracy("$(pkgroot)/data/F.csv", "grn_example2/Q_effects.tsv.gz"), digits = 2))
            Lasso+                                          $(round(GRN.accuracy("$(pkgroot)/data/F.csv", "grn_example3/Q_effects.tsv.gz"), digits = 2))
            Ridge                                           $(round(GRN.accuracy("$(pkgroot)/data/F.csv", "grn_example4/Q_effects.tsv.gz"), digits = 2))
            Spike and slab (no priors, plain)               $(round(GRN.accuracy("$(pkgroot)/data/F.csv", "grn_example5/Q_effects.tsv"), digits = 2))
            Spike and slab (no priors, binary)              $(round(GRN.accuracy("$(pkgroot)/data/F.csv", "grn_example6/Q_effects.bin"), digits = 2))
            Spike and slab (no priors, compr. binary)       $(round(GRN.accuracy("$(pkgroot)/data/F.csv", "grn_example7/Q_effects.bin.gz"), digits = 2))
            Spike and slab (continued from example 1)       $(round(GRN.accuracy("$(pkgroot)/data/F.csv", "grn_example8/Q_effects.tsv.gz"), digits = 2))
            "
