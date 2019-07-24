using Test, TruncatedNormal, SpecialFunctions

@test tnmom1c(0, -Inf, Inf) == 0
@test tnmom1c(0, -Inf, 0) ≈ tnmom1(-Inf, 0)
@test tnmom1c(0, 0, +Inf) ≈ tnmom1(0, +Inf)

@test tnmom1c(0, +Inf, +Inf) == +Inf
@test tnmom1c(0, -Inf, -Inf) == -Inf
@test tnmom1c(1, +Inf, +Inf) == +Inf
@test tnmom1c(1, -Inf, -Inf) == -Inf
@test tnmom1c(-Inf, +Inf, +Inf) == +Inf
@test tnmom1c(+Inf, -Inf, -Inf) == -Inf

@test isnan(tnmom1c(+Inf, +Inf, +Inf))
@test isnan(tnmom1c(-Inf, -Inf, -Inf))

@test isnan(tnmom1c(0, 1, 0))
@test isnan(tnmom1c(0, +Inf, 0))
@test isnan(tnmom1c(0, 0, -Inf))
@test isnan(tnmom1c(0, +Inf, -Inf))
@test isnan(tnmom1c(1, 1, 0))
@test isnan(tnmom1c(1, +Inf, 0))
@test isnan(tnmom1c(1, 0, -Inf))

for c = (0, 1), x = (0, NaN), y = (0, NaN)
    if isnan(x) || isnan(y)
        @test isnan(tnmom1c(c, x, y))
    end
end

for x = exp.(-100:100), c = -100:100
    @test tnmom1c(c, -x, x) ≈ -c
end

for x = -10:10, c = -10:10
    @test tnmom1c(c, x, +Inf) ≈ tnmom1(x, +Inf) - c
    @test tnmom1c(c, -Inf, x) ≈ tnmom1(-Inf, x) - c
end

#= tnmom1c(c, a, b) = 0 when c == tnmom1(a,b) =#
for a = exp.(-10:100), b = exp.(-10:100)
    a ≤ b || continue
    c = tnmom1(a, b)
    @test abs(tnmom1c(c, a, b)) ≤ max(1e-10, c * 1e-10)
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
@test tnmom1c(1, -1e200, 1e200) ≈ -1
@test tnmom1c(1e200, -1e200, 1e200) ≈ -1e200
@test tnmom1c(0, -1e6, -999000) ≈ tnmom1(-1e6, -999000)
@test tnmom1c(0, +1e6, +Inf) ≈ tnmom1(+1e6, +Inf)
@test tnmom1c(0, -Inf, -1e6) ≈ tnmom1(-Inf, -1e6)

@test tnmom1c(999000, 999000, 1e6) ≈ 1.00100100099899498898098100910e-6
@test tnmom1c(999500, 999000, 1e6) ≈ -499.999998998998999001005011019
@test tnmom1c(1e6, 999000, 1e6) ≈ -999.999998998998999001005011019
