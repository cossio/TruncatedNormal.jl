using Test, TruncatedNormal, SpecialFunctions

@test tnmom2c(0, -Inf, Inf) == 1
@test tnmom2c(0, -Inf, 0) == 1
@test tnmom2c(0, 0, +Inf) == 1

@test tnmom2c(0, +Inf, +Inf) == Inf
@test tnmom2c(0, -Inf, -Inf) == Inf

@test isnan(tnmom2c(0, 1, 0))
@test isnan(tnmom2c(0, +Inf, 0))
@test isnan(tnmom2c(0, 0, -Inf))

@test isnan(tnmom2c(0, NaN, NaN))
@test isnan(tnmom2c(0, 0, NaN))
@test isnan(tnmom2c(0, NaN, 0))

@test tnmom2c(10, 0, 10) ≈ 85.0423087839426928824034533318
@test tnmom2c(1, -3, 3) ≈ 1.97333692466254147658812248690
@test tnmom2c(1, 3, 4) ≈ 5.15893137101968604368598313652
@test tnmom2c(-1, 3, 4) ≈ 18.2007485133797706592985387254
@test tnmom2c(2, 3, 4) ≈ 1.63802279983964373587970534208
@test tnmom2c(-2, 3, 4) ≈ 27.7216570845598129671048165198

for x = -10:10
    @test tnmom2c(0, x, +Inf) ≈ 1 + x * √(2/π) / erfcx(+x/√2)
    @test tnmom2c(0, -Inf, x) ≈ 1 - x * √(2/π) / erfcx(-x/√2)
end

@test tnmom2c(0, -Inf, Inf, 0, 1) == 1

for x = -10:10
    @test tnmom2c(0, x, +Inf, 0, 1) ≈ tnmom2(x, +Inf)
    @test tnmom2c(0, -Inf, x, 0, 1) ≈ tnmom2(-Inf, x)
end

@test tnmom2c(0, -1e200, 1e200) ≈ 1
@test tnmom2c(0, -1e6, -999000) ≈ 9.9800100000199999999999799599399200002e11
@test tnmom2c(0, +1e6, +Inf) ≈ 1.00000000000199999999999800000e12
@test tnmom2c(0, -Inf, -1e6) ≈ 1.00000000000199999999999800000e12

@test tnmom2c(4, -10, 20, 3, 3//2) ≈ 3.24999999999999956053727803062
