using Test, Random, Zygote, Distributions, SpecialFunctions, FiniteDifferences
import TruncatedNormal as TN

for a = -10:10
    d = truncated(Normal(), a, Inf)
    @test TN.tnmean(a) ≈ mean(d)
    @test TN.tnstd(a)  ≈ std(d)
    @test TN.tnvar(a)  ≈ var(d)
    @test a ≤ TN.tnmean(a) < Inf
    @test 0 ≤ TN.tnvar(a) ≤ 1
end
@test 1e20 ≤ TN.tnmean(1e20) < Inf
@test_broken 0 ≤ TN.tnvar(1e80) ≤ 1

for a = -2:2
    @test mean(TN.randnt'(a) for _ = 1:10^6) ≈ TN.tnmean'(a) rtol=0.1
    z, da = TN.∇randnt(Random.GLOBAL_RNG, a)
    @test da ≈ erfc(z/√2) / erfc(a/√2) * exp((z^2 - a^2)/2)
end

@testset "randnt_half gradient" begin
    for μ = -1:1, σ = 1:2
        dμ_ = mean(gradient(TN.randnt_half, μ, σ)[1] for _ = 1:10^6)
        dσ_ = mean(gradient(TN.randnt_half, μ, σ)[2] for _ = 1:10^6)
        dμ, dσ = gradient(μ, σ) do μ, σ
            μ + σ * TN.tnmean(-μ/σ)
        end
        @test dμ ≈ dμ_ rtol=0.1
        @test dσ ≈ dσ_ rtol=0.1
    end
end

@testset "sqrt1half" begin
    @test TN.sqrt1half(5) ≈ 5.1925824035672520156
    @test TN.sqrt1half(0) == 1
    @test TN.sqrt1half(-1) == TN.sqrt1half(1) ≈ 1.6180339887498948482
    @test isnan(TN.sqrt1half(NaN))
    @test TN.sqrt1half(Inf) == TN.sqrt1half(-Inf) == Inf
    @test TN.sqrt1half(1e300) ≈ 1e300
end

@testset "randnt" begin
    @test TN.randnt(0) > 0
    @test TN.randnt(1e300) == 1e300
    @test TN.randnt(Inf) == Inf
    @test isnan(TN.randnt(NaN))
    @test TN.randnt(floatmax(Float64)) == floatmax(Float64)
end
