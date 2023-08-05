using Test: @test, @testset, @inferred
using TruncatedNormal: tnlogpdf
using Distributions: truncated, Normal, logpdf

@test tnlogpdf(-1, 0, 1) == -Inf
@test tnlogpdf(10, 0, 1) == -Inf
@test tnlogpdf(0.4, 0, 1) ≈ logpdf(truncated(Normal(), 0, 1), 0.4)

@inferred tnlogpdf(-1, 0, 1)
@inferred tnlogpdf(0.4, 0, 1)
