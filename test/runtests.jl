using GRN
using Random
using LinearAlgebra
using Statistics
using Test
using GZip
using DelimitedFiles
using Printf

# Helper: read an output file as a matrix with a given number of values per line.
# Handles both plain .tsv and gzip-compressed .tsv.gz files.
function readoutput(p, ncols)
    isfile(p) || error("Missing output file $p")
    io = endswith(p, ".gz") ? GZip.open(p) : open(p)
    try
        v = [parse(Float64, t) for ln in eachline(io) for t in split(ln, '\t') if !isempty(t)]
        return reshape(v, ncols, :)'
    finally
        close(io)
    end
end

@testset "GRN.jl" begin

    @testset "Input checks" begin
        out = joinpath(mktempdir(), "grn")
        @test_throws AssertionError grn(randn(5, 3), randn(6, 3), 2:2:4; outpath=out)
        @test_throws AssertionError grn(randn(10, 2), randn(10, 3), 4:2:2; outpath=out)
        @test_throws AssertionError grn(randn(10, 2), randn(10, 3), 0:2:4; outpath=out)
    end

    @testset "Prior selection" begin
        out = joinpath(mktempdir(), "grn")
        Y = randn(20, 3)
        @test_throws ArgumentError grn(Y, 2:2:4; priors=(Gprior=0.1,), outpath=out)
        @test_throws ArgumentError grn(Y, 2:2:4; priors=(Qprior="spikeslab",), outpath=out)
    end

    @testset "Scalar samplers" begin
        Random.seed!(1)
        # Gamma
        x = [GRN._rand_gamma(4.5) for _ in 1:100000]
        @test isapprox(mean(x), 4.5; rtol=0.02)
        @test isapprox(var(x), 4.5; rtol=0.05)
        # Chisq
        x = [GRN._rand_chisq(7.0) for _ in 1:100000]
        @test isapprox(mean(x), 7.0; rtol=0.02)
        @test isapprox(var(x), 14.0; rtol=0.05)
        # Beta
        x = [GRN._rand_beta(2.0, 5.0) for _ in 1:100000]
        @test isapprox(mean(x), 2.0 / 7.0; rtol=0.02)
        # Inverse Gaussian
        x = [GRN._rand_invgauss(2.0, 3.0) for _ in 1:100000]
        @test isapprox(mean(x), 2.0; rtol=0.05)
    end

    @testset "Simulate" begin
        # Default parameters
        sim = GRN.simulate(seed=1)
        @test size(sim.Y) == (200, 10)
        @test size(sim.M) == (200, 0)
        @test size(sim.G) == (10, 10)
        @test size(sim.Q) == (0, 10)
        @test length(sim.mu) == 10
        @test all(isfinite, sim.Y)

        # G is strictly upper triangular (G[i,j] = 0 for i >= j)
        lower_zero = all(sim.G[i, j] == 0.0 for j in 1:10 for i in j:10)
        @test lower_zero

        # With markers
        sim = GRN.simulate(n=100, m=5, u=20, seed=42)
        @test size(sim.Y) == (100, 5)
        @test size(sim.M) == (100, 20)
        @test size(sim.G) == (5, 5)
        @test size(sim.Q) == (20, 5)
        @test all(isfinite, sim.Y)
        @test all(isfinite, sim.M)

        # Reproducibility: same seed gives same data
        sim1 = GRN.simulate(n=50, m=4, u=8, seed=7)
        sim2 = GRN.simulate(n=50, m=4, u=8, seed=7)
        @test sim1.Y == sim2.Y
        @test sim1.M == sim2.M
        @test sim1.G == sim2.G

        # Different seeds give different data
        sim3 = GRN.simulate(n=50, m=4, u=8, seed=8)
        @test sim1.Y != sim3.Y

        # Model check: residuals Y - mu - Y*G - M*Q have sd close to noise_sd
        sim = GRN.simulate(n=500, m=5, u=10, seed=1, noise_sd=0.5)
        resid = sim.Y .- sim.mu' .- sim.Y * sim.G .- sim.M * sim.Q
        @test isapprox(std(resid), 0.5; rtol=0.15)

        # User-provided G (must be strictly upper triangular)
        G = zeros(3, 3)
        G[1, 2] = 0.5
        G[1, 3] = -0.3
        G[2, 3] = 0.7
        sim = GRN.simulate(n=50, m=3, G=G, seed=1)
        @test sim.G == G
        @test size(sim.Y) == (50, 3)

        # Non-upper-triangular G throws
        Gbad = zeros(3, 3)
        Gbad[2, 1] = 0.5
        @test_throws ArgumentError GRN.simulate(n=50, m=3, G=Gbad)

        # G with wrong dimensions throws
        @test_throws ArgumentError GRN.simulate(n=50, m=3, G=zeros(4, 4))

        # User-provided M: u inferred from size
        M = randn(80, 6)
        sim = GRN.simulate(n=80, M=M, m=4, seed=1)
        @test size(sim.M) == (80, 6)
        @test size(sim.Q) == (6, 4)

        # M with wrong number of rows throws
        @test_throws ArgumentError GRN.simulate(n=50, M=randn(60, 6), m=4)

        # User-provided mu
        sim = GRN.simulate(n=50, m=3, mu=[1.0, 2.0, 3.0], seed=1)
        @test sim.mu == [1.0, 2.0, 3.0]

        # mu with wrong length throws
        @test_throws ArgumentError GRN.simulate(n=50, m=3, mu=[1.0, 2.0])

        # Simulated data works with grn
        sim = GRN.simulate(n=50, m=3, u=4, seed=1)
        out = joinpath(mktempdir(), "grn")
        grn(sim.Y, sim.M, 10:10:50, outpath=out)
        @test isfile(joinpath(out, "G_effects.tsv.gz"))

        # Default: no fixed effects
        @test size(sim.X) == (50, 0)
        @test size(sim.Beta) == (0, 3)

        # With fixed effects
        X = [ones(50) randn(50)]
        sim = GRN.simulate(n=50, m=3, X=X, seed=1)
        @test size(sim.X) == (50, 2)
        @test size(sim.Beta) == (2, 3)
        @test all(sim.mu .== 0.0) # mu replaced by X * Beta

        # Model check with fixed effects
        sim = GRN.simulate(n=500, m=5, u=10, X=[ones(500) randn(500)], seed=1, noise_sd=0.5)
        resid = sim.Y .- sim.X * sim.Beta .- sim.Y * sim.G .- sim.M * sim.Q
        @test isapprox(std(resid), 0.5; rtol=0.15)

        # X with wrong rows throws
        @test_throws ArgumentError GRN.simulate(n=50, m=3, X=randn(60, 2))

        # Beta with wrong dims throws
        @test_throws ArgumentError GRN.simulate(n=50, m=3, X=randn(50, 2), Beta=zeros(3, 3))

        # User-provided Beta
        Beta = [1.0 2.0 3.0; 0.5 -0.5 0.2]
        sim = GRN.simulate(n=50, m=3, X=randn(50, 2), Beta=Beta, seed=1)
        @test sim.Beta == Beta
    end

    @testset "Fixed effects" begin
        # Input checks for X
        out = joinpath(mktempdir(), "grn")
        Y = randn(20, 3)
        @test_throws AssertionError grn(Y, 2:2:4; X=randn(10, 2), outpath=out)
        @test_throws AssertionError grn(Y, 2:2:4; X=[ones(20) fill(NaN, 20)], outpath=out)
        @test_throws AssertionError grn(Y, 2:2:4; X=zeros(20, 2), outpath=out)

        # Default X is ones(n, 1): intercept is captured as a fixed effect
        testdir = mktempdir()
        Random.seed!(10)
        sim = GRN.simulate(n=200, m=5, u=10, seed=1, effect_sd=0.8, noise_sd=0.25)
        grn(sim.Y, sim.M, 500:10:1000, outpath=joinpath(testdir, "default_x"))
        @test isfile(joinpath(testdir, "default_x", "base_fixedeffects.tsv.gz"))
        @test !isfile(joinpath(testdir, "default_x", "base_average.tsv.gz"))
        # 1 fixed effect (intercept) x 5 genes = 5 values per row, 51 iterations
        @test size(readoutput(joinpath(testdir, "default_x", "base_fixedeffects.tsv.gz"), 5)) == (51, 5)
        # The intercept captures the per-gene mean of the simulated data
        intercept_flat = vec(mean(readoutput(joinpath(testdir, "default_x", "base_fixedeffects.tsv.gz"), 5), dims=1))
        @test cor(intercept_flat, sim.mu) > 0.8

        # User-provided X with multiple fixed effects
        X = [ones(200) randn(200) rand(200) .> 0.5]
        sim = GRN.simulate(n=200, m=5, u=10, X=X, seed=1, effect_sd=0.8, noise_sd=0.25)
        grn(sim.Y, sim.M, 500:10:1000; X=sim.X, outpath=joinpath(testdir, "fe"))
        @test isfile(joinpath(testdir, "fe", "base_fixedeffects.tsv.gz"))
        @test !isfile(joinpath(testdir, "fe", "base_average.tsv.gz"))
        @test isfile(joinpath(testdir, "fe", "G_effects.tsv.gz"))
        @test isfile(joinpath(testdir, "fe", "Q_effects.tsv.gz"))
        # 3 fixed effects x 5 genes = 15 values per row, 51 iterations
        @test size(readoutput(joinpath(testdir, "fe", "base_fixedeffects.tsv.gz"), 15)) == (51, 15)

        # Recovery of fixed effects: posterior mean close to truth
        Beta_flat = vec(mean(readoutput(joinpath(testdir, "fe", "base_fixedeffects.tsv.gz"), 15), dims=1))
        Beta_hat = reshape(Beta_flat, 3, 5) # Beta_hat[k, c] = coefficient of fixed effect k for gene c
        @test cor(vec(Beta_hat), vec(sim.Beta)) > 0.9

        # Recovery of network with fixed effects present
        Gflat = vec(mean(readoutput(joinpath(testdir, "fe", "G_effects.tsv.gz"), 25), dims=1))
        Ghat = reshape(Gflat, 5, 5)' # Ghat[j, q] is the effect of gene q on gene j
        @test cor(vec(Ghat), vec(sim.G')) > 0.7

        # Fixed effects sampler is allocation-free
        sim2 = GRN.simulate(n=100, m=4, X=randn(100, 3), seed=1)
        d = GRN.create_databucket(sim2.Y, zeros(100, 0), (individuals=100, genes=4, loci=0),
            GRN.SpikeSlab(), GRN.SpikeSlab(), sim2.X)
        @test length(d.out.base.fixedeffects) == 4
        @test all(length(b) == 3 for b in d.out.base.fixedeffects)
        f_fe() = GRN._sample_fixedeffects!(d, 1)
        f_fe(); f_fe()
        @test (@allocated f_fe()) == 0

        # Without X: defaults to ones(n, 1) — intercept only
        d2 = GRN.create_databucket(sim2.Y, zeros(100, 0), (individuals=100, genes=4, loci=0),
            GRN.SpikeSlab(), GRN.SpikeSlab())
        @test length(d2.in.X) == 1
        @test all(d2.in.X[1] .== 1.0)
        @test length(d2.out.base.fixedeffects) == 4
        @test all(length(b) == 1 for b in d2.out.base.fixedeffects)

        # grn works with X and no genotypes (convenience method)
        sim3 = GRN.simulate(n=100, m=3, X=[ones(100) randn(100)], seed=2)
        out3 = joinpath(mktempdir(), "grn")
        grn(sim3.Y, 10:10:50; X=sim3.X, outpath=out3)
        @test isfile(joinpath(out3, "base_fixedeffects.tsv.gz"))
        @test !isfile(joinpath(out3, "base_average.tsv.gz"))
        @test !isfile(joinpath(out3, "Q_effects.tsv.gz"))
    end

    @testset "Priors are inserted" begin
        sim = GRN.simulate(n=50, m=3, u=4, seed=1)
        Y, M = sim.Y, sim.M
        d = GRN.create_databucket(Y, M, (individuals=50, genes=3, loci=4),
            GRN.SpikeSlab(), GRN.SpikeSlab())
        GRN.insert_priors!(d, (Rvar=0.5, Gpi=0.2, Gvar=0.3,
            Gvardf=6.0, Qpi=0.4, Rvardf=5.0))
        @test all(d.out.base.residual_variance .== 0.5)
        @test all(d.out.G.pi .== 0.2)
        @test all(d.out.G.var .== 0.3)
        @test all(d.out.G.varscale .== 0.3 * (6.0 - 2.0))
        @test all(d.out.Q.pi .== 0.4)
        @test all(d.out.base.residual_variancescale .== 0.5 * (5.0 - 2.0))
        # The G prior excludes the diagonal, the Q prior includes all loci
        @test [length(c) for c in d.out.G.candidates] == [2, 2, 2]
        @test [length(c) for c in d.out.Q.candidates] == [4, 4, 4]
    end

    @testset "Analysis and output" begin
        testdir = mktempdir()
        Random.seed!(2)
        sim = GRN.simulate(n=200, m=5, u=10, seed=1, effect_sd=0.8, noise_sd=0.25)
        Y, M, B = sim.Y, sim.M, sim.G

        # Run all prior combinations (compressed by default)
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "ss"))
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "lr"),
            priors=(Gprior=Laplace(), Qprior=Ridge()))
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "rl"),
            priors=(Gprior=Ridge(), Qprior=Laplace()))

        # Spike-and-slab output (compressed .tsv.gz, no indicator file)
        @test size(readoutput(joinpath(testdir, "ss", "G_effects.tsv.gz"), 25)) == (51, 25)
        @test !isfile(joinpath(testdir, "ss", "G_indicator.tsv.gz"))
        @test !isfile(joinpath(testdir, "ss", "G_effects.tsv"))
        @test size(readoutput(joinpath(testdir, "ss", "G_pi.tsv.gz"), 5)) == (51, 5)
        @test size(readoutput(joinpath(testdir, "ss", "G_variance.tsv.gz"), 5)) == (51, 5)
        @test size(readoutput(joinpath(testdir, "ss", "Q_effects.tsv.gz"), 50)) == (51, 50)
        @test size(readoutput(joinpath(testdir, "ss", "base_fixedeffects.tsv.gz"), 5)) == (51, 5)
        @test all(isfinite, readoutput(joinpath(testdir, "ss", "G_variance.tsv.gz"), 5))

        # Laplace output
        @test size(readoutput(joinpath(testdir, "lr", "G_effects.tsv.gz"), 25)) == (51, 25)
        @test size(readoutput(joinpath(testdir, "lr", "G_tau2.tsv.gz"), 25)) == (51, 25)
        @test size(readoutput(joinpath(testdir, "lr", "G_lambda2.tsv.gz"), 5)) == (51, 5)
        @test all(readoutput(joinpath(testdir, "lr", "G_tau2.tsv.gz"), 25) .> 0)
        @test all(readoutput(joinpath(testdir, "lr", "G_lambda2.tsv.gz"), 5) .> 0)

        # Ridge output
        @test size(readoutput(joinpath(testdir, "rl", "G_effects.tsv.gz"), 25)) == (51, 25)
        @test size(readoutput(joinpath(testdir, "rl", "G_variance.tsv.gz"), 5)) == (51, 5)

        # Analysis without genotypes writes no marker-gene files
        grn(Y, 500:10:1000, outpath=joinpath(testdir, "nom"))
        @test isfile(joinpath(testdir, "nom", "G_effects.tsv.gz"))
        @test !isfile(joinpath(testdir, "nom", "Q_effects.tsv.gz"))
        @test size(readoutput(joinpath(testdir, "nom", "G_effects.tsv.gz"), 25)) == (51, 25)

        # The network is recovered with the spike-and-slab prior
        Random.seed!(3)
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "rec"),
            priors=(Gpi=0.2, Qpi=0.05))
        Gflat = vec(mean(readoutput(joinpath(testdir, "rec", "G_effects.tsv.gz"), 25), dims=1))
        Ghat = reshape(Gflat, 5, 5)' # Ghat[j, q] is the effect of gene q on gene j
        @test cor(vec(Ghat), vec(B')) > 0.7
    end

    @testset "Compressed output (default)" begin
        testdir = mktempdir()
        Random.seed!(5)
        sim = GRN.simulate(n=200, m=5, u=10, seed=1)
        Y, M = sim.Y, sim.M

        # Run without specifying compress — default is compressed
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "default"))

        # .tsv.gz files exist, plain .tsv files do not
        @test isfile(joinpath(testdir, "default", "G_effects.tsv.gz"))
        @test isfile(joinpath(testdir, "default", "G_pi.tsv.gz"))
        @test isfile(joinpath(testdir, "default", "G_variance.tsv.gz"))
        @test isfile(joinpath(testdir, "default", "Q_effects.tsv.gz"))
        @test isfile(joinpath(testdir, "default", "base_fixedeffects.tsv.gz"))
        @test isfile(joinpath(testdir, "default", "base_residual_variance.tsv.gz"))
        @test !isfile(joinpath(testdir, "default", "G_effects.tsv"))
        @test !isfile(joinpath(testdir, "default", "Q_effects.tsv"))

        # Indicator files are never written (PIP derived from effects)
        @test !isfile(joinpath(testdir, "default", "G_indicator.tsv.gz"))
        @test !isfile(joinpath(testdir, "default", "Q_indicator.tsv.gz"))

        # Compressed files can be read back with correct dimensions
        @test size(readoutput(joinpath(testdir, "default", "G_effects.tsv.gz"), 25)) == (51, 25)
        @test size(readoutput(joinpath(testdir, "default", "G_pi.tsv.gz"), 5)) == (51, 5)
        @test size(readoutput(joinpath(testdir, "default", "G_variance.tsv.gz"), 5)) == (51, 5)
        @test size(readoutput(joinpath(testdir, "default", "Q_effects.tsv.gz"), 50)) == (51, 50)
        @test size(readoutput(joinpath(testdir, "default", "base_fixedeffects.tsv.gz"), 5)) == (51, 5)

        # All values are finite
        @test all(isfinite, readoutput(joinpath(testdir, "default", "G_effects.tsv.gz"), 25))
        @test all(isfinite, readoutput(joinpath(testdir, "default", "G_variance.tsv.gz"), 5))

        # posterior_mean works on .gz files
        pm = GRN.posterior_mean(joinpath(testdir, "default", "G_effects.tsv.gz"))
        @test length(pm) == 25
        @test all(isfinite, pm)

        # Compressed files are smaller than uncompressed
        Random.seed!(5)
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "plain"), compress=false)
        gz_size = filesize(joinpath(testdir, "default", "G_effects.tsv.gz"))
        plain_size = filesize(joinpath(testdir, "plain", "G_effects.tsv"))
        @test gz_size < plain_size
    end

    @testset "Plain output (compress=false)" begin
        testdir = mktempdir()
        Random.seed!(6)
        sim = GRN.simulate(n=200, m=5, u=10, seed=1)
        Y, M = sim.Y, sim.M

        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "plain"), compress=false)

        # Plain .tsv files exist, .tsv.gz files do not
        @test isfile(joinpath(testdir, "plain", "G_effects.tsv"))
        @test isfile(joinpath(testdir, "plain", "G_pi.tsv"))
        @test isfile(joinpath(testdir, "plain", "Q_effects.tsv"))
        @test isfile(joinpath(testdir, "plain", "base_fixedeffects.tsv"))
        @test !isfile(joinpath(testdir, "plain", "G_effects.tsv.gz"))
        @test !isfile(joinpath(testdir, "plain", "Q_effects.tsv.gz"))

        # Indicator files are never written
        @test !isfile(joinpath(testdir, "plain", "G_indicator.tsv"))

        # Plain files can be read back
        @test size(readoutput(joinpath(testdir, "plain", "G_effects.tsv"), 25)) == (51, 25)
        @test size(readoutput(joinpath(testdir, "plain", "G_pi.tsv"), 5)) == (51, 5)
        @test all(isfinite, readoutput(joinpath(testdir, "plain", "G_effects.tsv"), 25))

        # posterior_mean works on plain files
        pm = GRN.posterior_mean(joinpath(testdir, "plain", "G_effects.tsv"))
        @test length(pm) == 25
        @test all(isfinite, pm)
    end

    @testset "PIP and derive_indicator" begin
        testdir = mktempdir()
        Random.seed!(8)
        sim = GRN.simulate(n=200, m=5, u=10, seed=1)
        Y, M = sim.Y, sim.M

        # Run analysis (compressed by default, indicator files are never written)
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "run"))

        effects_path = joinpath(testdir, "run", "G_effects.tsv.gz")

        # --- posterior_inclusion_probability ---

        pip = GRN.posterior_inclusion_probability(effects_path)
        @test pip isa Vector{Float64}
        @test length(pip) == 25  # 5 genes x 5 predictors
        @test all(0.0 .<= pip .<= 1.0)

        # Verify PIP by manual computation: proportion of non-zero rows per column
        effects_mat = readoutput(effects_path, 25)  # 51 rows x 25 cols
        manual_pip = vec(mean(effects_mat .!= 0.0, dims=1))
        @test pip ≈ manual_pip

        # Diagonal always zero (G prior excludes diagonal)
        diag_idx = [1, 7, 13, 19, 25]
        @test all(pip[diag_idx] .== 0.0)

        # --- derive_indicator (PIP-based) ---

        # Default threshold 0.5
        ind = GRN.derive_indicator(effects_path)
        @test ind isa BitVector
        @test length(ind) == 25
        @test all(ind .== (pip .> 0.5))

        # threshold=0.0: included in at least one iteration
        ind_any = GRN.derive_indicator(effects_path; threshold=0.0)
        @test all(ind_any .== (pip .> 0.0))

        # threshold=0.8: high-confidence inclusions
        ind_strong = GRN.derive_indicator(effects_path; threshold=0.8)
        @test all(ind_strong .== (pip .> 0.8))

        # Diagonal always false
        @test all(.!ind[diag_idx])

        # Works on uncompressed files
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "plain"), compress=false)
        pip_plain = GRN.posterior_inclusion_probability(joinpath(testdir, "plain", "G_effects.tsv"))
        @test length(pip_plain) == 25
        ind_plain = GRN.derive_indicator(joinpath(testdir, "plain", "G_effects.tsv"))
        @test ind_plain isa BitVector
        @test all(ind_plain .== (pip_plain .> 0.5))

        # Marker-gene effects
        q_pip = GRN.posterior_inclusion_probability(joinpath(testdir, "run", "Q_effects.tsv.gz"))
        @test length(q_pip) == 50  # 5 genes x 10 markers
        @test all(0.0 .<= q_pip .<= 1.0)
    end

    @testset "Binary output (format=:bin)" begin
        testdir = mktempdir()
        Random.seed!(9)
        sim = GRN.simulate(n=200, m=5, u=10, seed=1, effect_sd=0.8, noise_sd=0.25)
        Y, M = sim.Y, sim.M

        # --- Compressed binary (.bin.gz) ---
        Random.seed!(9)
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "bin"),
            format=:bin, compress=true)

        # Correct extensions
        @test isfile(joinpath(testdir, "bin", "G_effects.bin.gz"))
        @test isfile(joinpath(testdir, "bin", "Q_effects.bin.gz"))
        @test isfile(joinpath(testdir, "bin", "G_pi.bin.gz"))
        @test isfile(joinpath(testdir, "bin", "base_fixedeffects.bin.gz"))
        @test !isfile(joinpath(testdir, "bin", "G_effects.tsv.gz"))
        @test !isfile(joinpath(testdir, "bin", "G_effects.bin"))

        # --- Uncompressed binary (.bin) ---
        Random.seed!(9)
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "bin_raw"),
            format=:bin, compress=false)

        @test isfile(joinpath(testdir, "bin_raw", "G_effects.bin"))
        @test isfile(joinpath(testdir, "bin_raw", "Q_effects.bin"))
        @test !isfile(joinpath(testdir, "bin_raw", "G_effects.bin.gz"))
        @test !isfile(joinpath(testdir, "bin_raw", "G_effects.tsv"))

        # --- read_binary: correct dimensions ---
        G_bin = GRN.read_binary(joinpath(testdir, "bin", "G_effects.bin.gz"))
        @test size(G_bin) == (51, 25)  # 51 iterations x 25 values (5x5)
        @test eltype(G_bin) == Float64
        @test all(isfinite, G_bin)

        Q_bin = GRN.read_binary(joinpath(testdir, "bin_raw", "Q_effects.bin"))
        @test size(Q_bin) == (51, 50)  # 51 iterations x 50 values (5x10)

        pi_bin = GRN.read_binary(joinpath(testdir, "bin", "G_pi.bin.gz"))
        @test size(pi_bin) == (51, 5)   # 51 iterations x 5 genes

        fe_bin = GRN.read_binary(joinpath(testdir, "bin", "base_fixedeffects.bin.gz"))
        @test size(fe_bin) == (51, 5)   # 51 iterations x 5 genes (1 fixed effect)

        # --- Binary vs TSV: same values ---
        Random.seed!(9)
        grn(Y, M, 500:10:1000, outpath=joinpath(testdir, "tsv"),
            format=:tsv, compress=true)

        G_tsv = readoutput(joinpath(testdir, "tsv", "G_effects.tsv.gz"), 25)
        @test size(G_tsv) == size(G_bin)
        @test G_tsv ≈ G_bin rtol=1e-9  # TSV has %.10g precision

        # --- posterior_mean works with binary ---
        pm_bin = GRN.posterior_mean(joinpath(testdir, "bin", "G_effects.bin.gz"))
        @test length(pm_bin) == 25
        @test all(isfinite, pm_bin)

        pm_bin_raw = GRN.posterior_mean(joinpath(testdir, "bin_raw", "G_effects.bin"))
        @test pm_bin ≈ pm_bin_raw  # compressed and uncompressed give same result

        # Compare with TSV posterior_mean
        pm_tsv = GRN.posterior_mean(joinpath(testdir, "tsv", "G_effects.tsv.gz"))
        @test pm_bin ≈ pm_tsv rtol=1e-9

        # --- posterior_inclusion_probability works with binary ---
        pip_bin = GRN.posterior_inclusion_probability(
            joinpath(testdir, "bin", "G_effects.bin.gz"))
        @test length(pip_bin) == 25
        @test all(0.0 .<= pip_bin .<= 1.0)

        # Diagonal always zero
        diag_idx = [1, 7, 13, 19, 25]
        @test all(pip_bin[diag_idx] .== 0.0)

        # Compare with TSV PIP
        pip_tsv = GRN.posterior_inclusion_probability(
            joinpath(testdir, "tsv", "G_effects.tsv.gz"))
        @test pip_bin ≈ pip_tsv

        # --- derive_indicator works with binary ---
        ind_bin = GRN.derive_indicator(joinpath(testdir, "bin", "G_effects.bin.gz"))
        @test ind_bin isa BitVector
        @test length(ind_bin) == 25
        @test all(ind_bin .== (pip_bin .> 0.5))

        # --- accuracy works with binary ---
        # B.csv is 10x10 but sim.G is 5x5; use a matching truth file
        B_truth = joinpath(testdir, "B_truth.csv")
        writedlm(B_truth, sim.G', ',')
        acc_bin = GRN.accuracy(B_truth, joinpath(testdir, "bin", "G_effects.bin.gz"))
        acc_tsv = GRN.accuracy(B_truth, joinpath(testdir, "tsv", "G_effects.tsv.gz"))
        @test acc_bin isa Float64
        @test acc_bin ≈ acc_tsv rtol=1e-9

        # --- Invalid format throws ---
        @test_throws ArgumentError grn(Y, M, 2:2:4,
            outpath=joinpath(testdir, "bad"), format=:csv)

        # --- Corrupt binary file throws ---
        corrupt = joinpath(testdir, "corrupt.bin")
        write(corrupt, "NOT A GRN BINARY FILE AT ALL")
        @test_throws ErrorException GRN.read_binary(corrupt)

        # --- Binary write allocates minimally ---
        # The binary write path should not allocate per element
        d = GRN.create_databucket(Y, M, (individuals=200, genes=5, loci=10),
            GRN.SpikeSlab(), GRN.SpikeSlab())
        ofiles = GRN.makeout(joinpath(testdir, "alloc_test"),
            GRN._outputs(d, 10); force=true, format=:bin, compress=false)
        f_write() = GRN.write_output(ofiles[1])  # G_effects
        f_write(); f_write()
        alloc = @allocated f_write()
        GRN.close_outputs(ofiles)
        # Binary write of 25 Float64 should allocate very little
        # (buffer flush may allocate when crossing _BUFSIZE threshold)
        @test alloc < 4096  # generous bound for buffer management
    end

    @testset "Preliminary results" begin
        testdir = mktempdir()
        Random.seed!(11)
        nr, nc = 40, 25
        V = randn(nr, nc)
        V[rand(nr, nc) .< 0.4] .= 0.0           # spike-and-slab-like zeros
        fmt(r) = join((@sprintf("%.10g", x) for x in V[r, :]), '\t')
        ref(k) = (vec(mean(V[1:k, :], dims=1)), vec(mean(V[1:k, :] .!= 0.0, dims=1)))

        # Finished plain file: equals posterior_mean and PIP
        p = joinpath(testdir, "full.tsv")
        write(p, join((fmt(r) * "\n" for r in 1:nr)))
        res = preliminary_posterior(p)
        @test res.iterations == nr
        @test res.mean ≈ GRN.posterior_mean(p)
        @test res.pip ≈ GRN.posterior_inclusion_probability(p)

        # Last row cut mid-number, no end-of-line: the row is dropped
        p = joinpath(testdir, "cut.tsv")
        lastrow = fmt(nr)
        write(p, join((fmt(r) * "\n" for r in 1:nr-1)) * lastrow[1:div(length(lastrow), 2)])
        res = preliminary_posterior(p)
        @test res.iterations == nr - 1
        @test res.mean ≈ ref(nr - 1)[1]
        @test res.pip ≈ ref(nr - 1)[2]

        # Complete last row without end-of-line is also dropped (cannot be
        # told apart from a row cut exactly at a value boundary)
        p = joinpath(testdir, "noeol.tsv")
        write(p, join((fmt(r) for r in 1:nr), '\n'))
        @test preliminary_posterior(p).iterations == nr - 1

        # CRLF line endings and blank lines
        p = joinpath(testdir, "crlf.tsv")
        write(p, "\r\n" * join((fmt(r) * "\r\n" for r in 1:nr)) * "\r\n")
        res = preliminary_posterior(p)
        @test res.iterations == nr
        @test res.mean ≈ ref(nr)[1]

        # Rows longer than the 1 MB read chunk: fields split across chunks
        W = randn(6, 60_000)
        p = joinpath(testdir, "wide.tsv")
        write(p, join((join((@sprintf("%.10g", x) for x in W[r, :]), '\t') * "\n" for r in 1:6)))
        @test filesize(p) > 2 * GRN._BUFSIZE
        res = preliminary_posterior(p)
        @test res.iterations == 6
        @test res.mean ≈ vec(mean(W, dims=1))

        # Gzip stream without trailer (truncated): complete rows before the cut
        p = joinpath(testdir, "full.tsv.gz")
        GZip.open(p, "w", 1) do io
            for r in 1:nr
                write(io, fmt(r) * "\n")
            end
        end
        @test preliminary_posterior(p).iterations == nr
        bytes = read(p)
        q = joinpath(testdir, "trunc.tsv.gz")
        write(q, bytes[1:end-30])
        res = preliminary_posterior(q)
        @test 0 < res.iterations < nr
        @test res.mean ≈ ref(res.iterations)[1]
        @test res.pip ≈ ref(res.iterations)[2]

        # Empty file and malformed values
        p = joinpath(testdir, "empty.tsv")
        write(p, "")
        res = preliminary_posterior(p)
        @test res.iterations == 0 && isempty(res.mean) && isempty(res.pip)
        p = joinpath(testdir, "bad.tsv")
        write(p, "1.0\tx\n")
        @test_throws ErrorException preliminary_posterior(p)
        p = joinpath(testdir, "ragged.tsv")
        write(p, "1.0\t2.0\n1.0\n")
        @test_throws ErrorException preliminary_posterior(p)

        # Binary: incomplete last row, missing header, truncated gzip
        for compress in (false, true)
            ext = compress ? ".bin.gz" : ".bin"
            out = joinpath(testdir, "bin$compress")
            vals = zeros(nc)
            files = GRN.makeout(out, (x=vals,); compress=compress, format=:bin)
            for r in 1:nr
                vals .= V[r, :]
                GRN.write_output(files[1])
            end
            GRN.close_outputs(files)
            p = joinpath(out, "x$ext")
            res = preliminary_posterior(p)
            @test res.iterations == nr
            @test res.mean ≈ ref(nr)[1]
            @test res.pip ≈ ref(nr)[2]

            if !compress
                bytes = read(p)
                q = joinpath(testdir, "partialrow.bin")
                write(q, bytes[1:end-3*8])
                @test preliminary_posterior(q).iterations == nr - 1
                write(q, bytes[1:10])
                @test preliminary_posterior(q).iterations == 0
            else
                bytes = read(p)
                q = joinpath(testdir, "trunc.bin.gz")
                write(q, bytes[1:end-30])
                res = preliminary_posterior(q)
                @test 0 < res.iterations <= nr
                @test res.mean ≈ ref(res.iterations)[1]
            end
        end

        # Files of a running analysis: written through makeout and still open
        for (compress, format) in ((true, :tsv), (false, :tsv), (true, :bin), (false, :bin))
            out = joinpath(testdir, "running_$(compress)_$format")
            vals = zeros(nc)
            files = GRN.makeout(out, (x=vals,); compress=compress, format=format)
            p = joinpath(out, only(readdir(out)))
            try
                @test preliminary_posterior(p).iterations == 0   # nothing flushed yet
                for r in 1:nr
                    vals .= V[r, :]
                    GRN.write_output(files[1])
                end
                # Rows reach the disk only after the 1 MB buffer and the
                # stream's own buffer (IOStream or zlib) are flushed
                GRN.flush_outputs(files)
                foreach(o -> flush(o.io), files)
                res = preliminary_posterior(p)
                @test res.iterations == nr
                @test res.mean ≈ ref(nr)[1]
            finally
                GRN.close_outputs(files)
            end
        end
    end

    @testset "Preliminary columns" begin
        testdir = mktempdir()
        Random.seed!(12)
        nr, nc = 30, 17
        V = randn(nr, nc)
        V[rand(nr, nc) .< 0.4] .= 0.0

        # Write V in all four formats through makeout, as grn does
        paths = String[]
        for (compress, format) in ((false, :tsv), (true, :tsv), (false, :bin), (true, :bin))
            out = joinpath(testdir, "cols_$(compress)_$format")
            vals = zeros(nc)
            files = GRN.makeout(out, (x=vals,); compress=compress, format=format)
            for r in 1:nr
                vals .= V[r, :]
                GRN.write_output(files[1])
            end
            GRN.close_outputs(files)
            push!(paths, joinpath(out, only(readdir(out))))
        end

        for p in paths
            # Single column, first and last column
            @test isapprox(preliminary_column(p, 5), V[:, 5]; rtol=1e-9)
            @test isapprox(preliminary_column(p, 1), V[:, 1]; rtol=1e-9)
            @test isapprox(preliminary_column(p, nc), V[:, nc]; rtol=1e-9)

            # Several columns in the given order, with duplicates
            idx = [9, 2, 9, 17]
            tr = preliminary_column(p, idx)
            @test size(tr) == (nr, 4)
            @test isapprox(tr, V[:, idx]; rtol=1e-9)

            # Out-of-range and non-positive columns
            @test_throws ArgumentError preliminary_column(p, nc + 1)
            @test_throws ArgumentError preliminary_column(p, 0)
        end

        # Unfinished files: only complete iterations
        tsv, bin = paths[1], paths[3]
        bytes = read(tsv)
        q = joinpath(testdir, "cut.tsv")
        write(q, bytes[1:end-5])                        # last row cut, no end-of-line
        @test isapprox(preliminary_column(q, 3), V[1:nr-1, 3]; rtol=1e-9)
        @test isapprox(preliminary_column(q, [nc, 1]), V[1:nr-1, [nc, 1]]; rtol=1e-9)

        bytes = read(bin)
        q = joinpath(testdir, "cut.bin")
        write(q, bytes[1:end-4*8])                      # last row missing its last 4 values
        @test isapprox(preliminary_column(q, 2), V[1:nr-1, 2]; rtol=1e-9)
        @test isapprox(preliminary_column(q, nc), V[1:nr-1, nc]; rtol=1e-9)
        write(q, bytes[1:10])                           # header incomplete
        @test preliminary_column(q, 2) == Float64[]
        @test size(preliminary_column(q, [1, 2])) == (0, 2)

        for p in (paths[2], paths[4])                   # gzip without trailer
            bytes = read(p)
            q = joinpath(testdir, "trunc_" * basename(p))
            write(q, bytes[1:end-30])
            tr = preliminary_column(q, [4, 11])
            k = size(tr, 1)
            @test 0 < k <= nr
            @test isapprox(tr, V[1:k, [4, 11]]; rtol=1e-9)
            @test k == preliminary_posterior(q).iterations
        end

        # Rows wider than the 1 MB read chunk
        W = randn(4, 60_000)
        q = joinpath(testdir, "wide.tsv")
        write(q, join((join((@sprintf("%.17g", x) for x in W[r, :]), '\t') * "\n" for r in 1:4)))
        @test preliminary_column(q, [1, 30_000, 60_000]) == W[:, [1, 30_000, 60_000]]

        # Column layout of an effects file from grn: effect of q on gene c
        # is column (c - 1) * m + q
        m = 5
        sim = GRN.simulate(n=50, m=m, u=0, seed=3)
        out = joinpath(testdir, "grnrun")
        redirect_stdout(devnull) do
            grn(sim.Y, 10:2:30; outpath=out, format=:bin, compress=false)
        end
        Gall = GRN.read_binary(joinpath(out, "G_effects.bin"))
        c, q = 4, 2
        @test preliminary_column(joinpath(out, "G_effects.bin"), (c - 1) * m + q) ==
            Gall[:, (c - 1) * m + q]
        @test all(iszero, preliminary_column(joinpath(out, "G_effects.bin"), (c - 1) * m + c))
    end

    @testset "Effective sample size" begin
        # Normal quantile function
        @test GRN._qnorm(0.5) == 0.0
        @test GRN._qnorm(0.975) ≈ 1.959963984540054 rtol = 1e-14
        @test GRN._qnorm(1e-10) ≈ -6.361340902404056 rtol = 1e-14
        @test GRN._qnorm(0.3) ≈ -GRN._qnorm(0.7) rtol = 1e-14

        # Reference values from MCMCDiagnosticTools.ess(...; maxlag=typemax(Int))
        # on a deterministic two-chain input with ties at zero
        t = 1:301
        x = hcat(sin.(0.3 .* t) .+ 0.5 .* cos.(37.0 .* t),
            cos.(0.25 .* t) .+ 0.5 .* sin.(41.0 .* t) .+ 0.1)
        x[x .< -0.6] .= 0.0
        @test ess(x) ≈ 144.54420517595005 rtol = 1e-10
        @test ess(x; kind=:tail) ≈ 372.493765673391 rtol = 1e-10
        @test ess(x; kind=:basic) ≈ 131.37705984571465 rtol = 1e-10

        # Independent draws: ESS close to the number of draws
        rng = Xoshiro(21)
        iid = randn(rng, 4000, 4)
        @test 0.85 * 16000 < ess(iid) < 1.15 * 16000

        # AR(1) with rho = 0.9: ESS close to S * (1 - rho) / (1 + rho)
        ar = zeros(5000, 4)
        for j in 1:4
            ar[1, j] = randn(rng)
            for i in 2:5000
                ar[i, j] = 0.9 * ar[i-1, j] + sqrt(1 - 0.81) * randn(rng)
            end
        end
        expected = 20000 * 0.1 / 1.9
        @test 0.75 * expected < ess(ar; kind=:basic) < 1.25 * expected

        # Bulk ESS is rank based: invariant to monotone transformations
        @test ess(exp.(ar)) ≈ ess(ar)

        # A drifting chain has a low ESS because the chain is split
        drift = collect(range(0, 10, length=2000)) .+ 0.1 .* randn(rng, 2000)
        @test ess(drift) < 50

        # Vector input, Bool input (inclusion indicators), undefined cases
        @test ess(ar[:, 1]) == ess(reshape(ar[:, 1], :, 1))
        @test isfinite(ess(ar[:, 1] .> 0.0; kind=:basic))
        @test isnan(ess(zeros(100)))                    # constant
        @test isnan(ess(randn(rng, 9)))                 # too short
        @test isnan(ess([1.0, 2.0, NaN, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]))
        @test_throws ArgumentError ess(iid; kind=:foo)

        # From files, one or several chains
        testdir = mktempdir()
        V = randn(rng, 300, 6)
        paths = String[]
        for (k, nr) in enumerate((300, 250))           # second chain is shorter
            out = joinpath(testdir, "chain$k")
            vals = zeros(6)
            files = GRN.makeout(out, (x=vals,); compress=true, format=:bin)
            for r in 1:nr
                vals .= V[r, :] .+ k
                GRN.write_output(files[1])
            end
            GRN.close_outputs(files)
            push!(paths, joinpath(out, "x.bin.gz"))
        end
        @test ess(paths[1], 3) == ess(V[:, 3] .+ 1)
        @test ess(paths[1], [3, 5]; kind=:tail) == [ess(V[:, c] .+ 1; kind=:tail) for c in (3, 5)]
        @test ess(paths, 2) == ess(hcat(V[1:250, 2] .+ 1, V[1:250, 2] .+ 2))
    end

    @testset "Allocations" begin
        # The per-gene Gibbs updates must not allocate
        Random.seed!(4)
        sim = GRN.simulate(n=100, m=4, u=6, seed=1)
        Y, M = sim.Y, sim.M
        d = GRN.create_databucket(Y, M, (individuals=100, genes=4, loci=6),
            GRN.SpikeSlab(), GRN.SpikeSlab())
        ssqs = (M=[dot(m, m) for m in d.in.M], Y=[dot(y, y) for y in d.in.Y])
        rv = d.out.base.residual_variance

        f() = begin
            for c in 1:4
                GRN._sample_fixedeffects!(d, c)
                GRN._sample_residualvariance!(d, c)
                GRN._sample_effects_gene!(d.out.G, c, d.in.Y, ssqs.Y, d.in.R, rv)
                GRN._sample_hyperparameters_gene!(d.out.G, c, rv)
                GRN._sample_effects_gene!(d.out.Q, c, d.in.M, ssqs.M, d.in.R, rv)
                GRN._sample_hyperparameters_gene!(d.out.Q, c, rv)
            end
        end
        f(); f()
        @test (@allocated f()) == 0

        # One iteration only allocates the tasks spawned by @threads
        f2() = GRN._sample_iteration!(d, ssqs, 100, 6)
        f2(); f2()
        @test (@allocated f2()) < 6 * 4096
    end
end
