using Test, TruncatedNormal, SpecialFunctions

@test tnmom2(-Inf, Inf) == 1
@test tnmom2(+Inf, +Inf) == Inf
@test tnmom2(-Inf, -Inf) == Inf
@test tnmom2(-Inf, 0) == 1
@test tnmom2(0, +Inf) == 1

@test isnan(tnmom2(1, 0))
@test isnan(tnmom2(+Inf, 0))
@test isnan(tnmom2(0, -Inf))
@test isnan(tnmom2(+Inf, -Inf))
@test isnan(tnmom2(NaN, NaN))
@test isnan(tnmom2(0, NaN))
@test isnan(tnmom2(NaN, 0))

for x = -10:10
    @test tnmom2(x, +Inf) ≈ 1 + x * √(2/π) / erfcx(+x/√2)
    @test tnmom2(-Inf, x) ≈ 1 - x * √(2/π) / erfcx(-x/√2)
end

@test tnmom2(-Inf, Inf, 0, 1) == 1

for x = -10:10
    @test tnmom2(x, +Inf, 0, 1) ≈ tnmom2(x, +Inf)
    @test tnmom2(-Inf, x, 0, 1) ≈ tnmom2(-Inf, x)
end

@test tnmom2(-1e200, 1e200) ≈ 1
@test tnmom2(-1e6, -999000) ≈ 9.9800100000199999999999799599399200002e11
@test tnmom2(+1e6, +Inf) ≈ 1.00000000000199999999999800000e12
@test tnmom2(-Inf, -1e6) ≈ 1.00000000000199999999999800000e12

@test tnmom2(-Inf, Inf, 2, 0) ≈ 4
@test tnmom2(0, 20, 3, 3//2) ≈ 11.4986153820554548159603543747
@test tnmom2(10, 20, 3, 3//2) ≈ 106.111726253622178129044932133
@test tnmom2(-10, 20, 3, 3//2) ≈ 11.2499999999999997949173964131
