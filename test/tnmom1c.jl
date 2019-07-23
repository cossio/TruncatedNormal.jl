using Test, TruncatedNormal, SpecialFunctions

@test tnmom1c(0, -Inf, Inf) == 0
@test tnmom1c(0, -Inf, 0) ≈ tnmom1(-Inf, 0)
@test tnmom1c(0, 0, +Inf) ≈ tnmom1(0, +Inf)

@test tnmom1c(0, +Inf, +Inf) == +Inf
@test tnmom1c(0, -Inf, -Inf) == -Inf

@test isnan(tnmom1c(0, 1, 0))
@test isnan(tnmom1c(0, +Inf, 0))
@test isnan(tnmom1c(0, 0, -Inf))

@test isnan(tnmom1c(0, NaN, NaN))
@test isnan(tnmom1c(0, 0, NaN))
@test isnan(tnmom1c(0, NaN, 0))

for x = exp.(-100:100), c = -100:100
    @test tnmom1c(c, -x, x) ≈ -c
end

for x = -10:10
    @test tnmom1c(0, x, +Inf) ≈ tnmom1(x, +Inf)
    @test tnmom1c(0, -Inf, x) ≈ tnmom1(-Inf, x)
end

@test tnmom1c(0, -Inf, Inf, 0, 1) ≈ 0

for x = -10:10
    @test tnmom1c(0, x, +Inf, 0, 1) ≈ tnmom1(x, +Inf)
    @test tnmom1c(0, -Inf, x, 0, 1) ≈ tnmom1(-Inf, x)
end

@test tnmom1c(10, 0, 10) ≈ -9.20211543919713464412024961258
@test tnmom1c(1, 3, 4) ≈ 2.26045428559002115390313889722
@test tnmom1c(-1, 3, 4) ≈ 4.26045428559002115390313889722
@test tnmom1c(2, 3, 4) ≈ 1.26045428559002115390313889722
@test tnmom1c(-2, 3, 4) ≈ 5.26045428559002115390313889722
@test tnmom1c(0, -1e200, 1e200) ≈ 0
@test tnmom1c(0, -1e6, -999000) ≈ tnmom1(-1e6, -999000)
@test tnmom1c(0, +1e6, +Inf) ≈ tnmom1(+1e6, +Inf)
@test tnmom1c(0, -Inf, -1e6) ≈ tnmom1(-Inf, -1e6)

@test tnmom1c(999000, 999000, 1e6) ≈ 1.00100100099899498898098100910e-6
