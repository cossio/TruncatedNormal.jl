using Test, SpecialFunctions
import TruncatedNormal as TN

@test TN.tnmom2c(0, -Inf, Inf) == 1
@test TN.tnmom2c(0, -Inf, 0) == 1
@test TN.tnmom2c(0, 0, +Inf) == 1

@test TN.tnmom2c(0, +Inf, +Inf) == Inf
@test TN.tnmom2c(0, -Inf, -Inf) == Inf

@test isnan(TN.tnmom2c(0, 1, 0))
@test isnan(TN.tnmom2c(0, +Inf, 0))
@test isnan(TN.tnmom2c(0, 0, -Inf))

@test isnan(TN.tnmom2c(0, NaN, NaN))
@test isnan(TN.tnmom2c(0, 0, NaN))
@test isnan(TN.tnmom2c(0, NaN, 0))

@test TN.tnmom2c(10, 0, 10) ≈ 85.0423087839426928824034533318
@test TN.tnmom2c(1, -3, 3) ≈ 1.97333692466254147658812248690
@test TN.tnmom2c(1, 3, 4) ≈ 5.15893137101968604368598313652
@test TN.tnmom2c(-1, 3, 4) ≈ 18.2007485133797706592985387254
@test TN.tnmom2c(2, 3, 4) ≈ 1.63802279983964373587970534208
@test TN.tnmom2c(-2, 3, 4) ≈ 27.7216570845598129671048165198

for x = -10:10
    @test TN.tnmom2c(0, x, +Inf) ≈ 1 + x * √(2/π) / erfcx(+x/√2)
    @test TN.tnmom2c(0, -Inf, x) ≈ 1 - x * √(2/π) / erfcx(-x/√2)
end

@test TN.tnmom2c(0, -Inf, Inf, 0, 1) == 1

for x = -10:10
    @test TN.tnmom2c(0, x, +Inf, 0, 1) ≈ TN.tnmom2(x, +Inf)
    @test TN.tnmom2c(0, -Inf, x, 0, 1) ≈ TN.tnmom2(-Inf, x)
end

@test TN.tnmom2c(0, -1e200, 1e200) ≈ 1
@test TN.tnmom2c(0, -1e6, -999000) ≈ 9.9800100000199999999999799599399200002e11
@test TN.tnmom2c(0, +1e6, +Inf) ≈ 1.00000000000199999999999800000e12
@test TN.tnmom2c(0, -Inf, -1e6) ≈ 1.00000000000199999999999800000e12

@test TN.tnmom2c(4, -10, 20, 3, 3//2) ≈ 3.24999999999999956053727803062
