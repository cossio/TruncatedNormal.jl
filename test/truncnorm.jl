using Test, Random, Zygote, Distributions, SpecialFunctions, FiniteDifferences
import TruncatedNormal as TruncNorm

for a = -10:10
    d = truncated(Normal(), a, Inf)
    @test TruncNorm.tnmean(a) ≈ mean(d)
    @test TruncNorm.tnstd(a)  ≈ std(d)
    @test TruncNorm.tnvar(a)  ≈ var(d)
    @test a ≤ TruncNorm.tnmean(a) < Inf
    @test 0 ≤ TruncNorm.tnvar(a) ≤ 1
end
@test 1e20 ≤ TruncNorm.tnmean(1e20) < Inf
@test_broken 0 ≤ TruncNorm.tnvar(1e80) ≤ 1

for a = -2:2
    @test mean(TruncNorm.randnt'(a) for _ = 1:10^6) ≈ TruncNorm.tnmean'(a) rtol=0.1
    z, da = TruncNorm.∇randnt(Random.GLOBAL_RNG, a)
    @test da ≈ erfc(z/√2) / erfc(a/√2) * exp((z^2 - a^2)/2)
end

@testset "randnt_half gradient" begin
    for μ = -1:1, σ = 1:2
        dμ_ = mean(gradient(TruncNorm.randnt_half, μ, σ)[1] for _ = 1:10^6)
        dσ_ = mean(gradient(TruncNorm.randnt_half, μ, σ)[2] for _ = 1:10^6)
        dμ, dσ = gradient(μ, σ) do μ, σ
            μ + σ * TruncNorm.tnmean(-μ/σ)
        end
        @test dμ ≈ dμ_ rtol=0.1
        @test dσ ≈ dσ_ rtol=0.1
    end
end

@testset "sqrt1half" begin
    @test TruncNorm.sqrt1half(5) ≈ 5.1925824035672520156
    @test TruncNorm.sqrt1half(0) == 1
    @test TruncNorm.sqrt1half(-1) == TruncNorm.sqrt1half(1) ≈ 1.6180339887498948482
    @test isnan(TruncNorm.sqrt1half(NaN))
    @test TruncNorm.sqrt1half(Inf) == TruncNorm.sqrt1half(-Inf) == Inf
    @test TruncNorm.sqrt1half(1e300) ≈ 1e300
end

@testset "randnt" begin
    @test TruncNorm.randnt(0) > 0
    @test TruncNorm.randnt(1e300) == 1e300
    @test TruncNorm.randnt(Inf) == Inf
    @test isnan(TruncNorm.randnt(NaN))
    @test TruncNorm.randnt(floatmax(Float64)) == floatmax(Float64)
end
