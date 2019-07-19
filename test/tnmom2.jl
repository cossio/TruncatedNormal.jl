using Test, TruncatedNormal, SpecialFunctions

@test tnmom2(-Inf, Inf) == 1
@test tnmom2(-Inf, 0) == 1
@test tnmom2(0, +Inf) == 1
@test tnmom2(+Inf, +Inf) == Inf
@test tnmom2(-Inf, -Inf) == Inf

@test isnan(tnmom2(1, 0))
@test isnan(tnmom2(+Inf, 0))
@test isnan(tnmom2(0, -Inf))

@test isnan(tnmom2(NaN, NaN))
@test isnan(tnmom2(0, NaN))
@test isnan(tnmom2(NaN, 0))

@test tnmom2(-Inf, Inf, 0, 1) == 1

@test tnmom2(0, +Inf) == 1
@test tnmom2(-Inf, 0) == 1

for x = -10:10
    @test tnmom2(x, +Inf) ≈ 1 + x * √(2/π) / erfcx(+x/√2)
    @test tnmom2(-Inf, x) ≈ 1 - x * √(2/π) / erfcx(-x/√2)
end

for x = -10:10
    @test tnmom2(x, +Inf) == tnmom2(x, +Inf, 0, 1)
    @test tnmom2(-Inf, x) == tnmom2(-Inf, x, 0, 1)
end

@test tnmom2(-1e200, 1e200) ≈ 1
@test tnmom2(-1e6, -999000) ≈ 9.9800100000199999999999799599399200002e11
