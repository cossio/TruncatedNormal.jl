using Test, SpecialFunctions
import TruncatedNormal as TN

@test TN.tnmom1c(0, -Inf, Inf) == 0
@test TN.tnmom1c(0, -Inf, 0) ≈ TN.tnmean(-Inf, 0)
@test TN.tnmom1c(0, 0, +Inf) ≈ TN.tnmean(0, +Inf)

@test TN.tnmom1c(0, +Inf, +Inf) == +Inf
@test TN.tnmom1c(0, -Inf, -Inf) == -Inf
@test TN.tnmom1c(1, +Inf, +Inf) == +Inf
@test TN.tnmom1c(1, -Inf, -Inf) == -Inf
@test TN.tnmom1c(-Inf, +Inf, +Inf) == +Inf
@test TN.tnmom1c(+Inf, -Inf, -Inf) == -Inf

@test isnan(TN.tnmom1c(+Inf, +Inf, +Inf))
@test isnan(TN.tnmom1c(-Inf, -Inf, -Inf))

@test isnan(TN.tnmom1c(0, 1, 0))
@test isnan(TN.tnmom1c(0, +Inf, 0))
@test isnan(TN.tnmom1c(0, 0, -Inf))
@test isnan(TN.tnmom1c(0, +Inf, -Inf))
@test isnan(TN.tnmom1c(1, 1, 0))
@test isnan(TN.tnmom1c(1, +Inf, 0))
@test isnan(TN.tnmom1c(1, 0, -Inf))

for c = (0, 1), x = (0, NaN), y = (0, NaN)
    if isnan(x) || isnan(y)
        @test isnan(TN.tnmom1c(c, x, y))
    end
end

for x = exp.(-100:100), c = -100:100
    @test TN.tnmom1c(c, -x, x) ≈ -c
end

for x = -10:10, c = -10:10
    @test TN.tnmom1c(c, x, +Inf) ≈ TN.tnmean(x, +Inf) - c
    @test TN.tnmom1c(c, -Inf, x) ≈ TN.tnmean(-Inf, x) - c
end

#= TN.tnmom1c(c, a, b) = 0 when c == tnmean(a,b) =#
for a = exp.(-10:100), b = exp.(-10:100)
    a ≤ b || continue
    c = TN.tnmean(a, b)
    @test abs(TN.tnmom1c(c, a, b)) ≤ max(1e-10, c * 1e-10)
    @test TN.tnmom1c(a, a, b) ≥ 0
    @test TN.tnmom1c(b, a, b) ≤ 0
end

@test TN.tnmom1c(0, -Inf, Inf, 0, 1) ≈ 0

for x = -10:10
    @test TN.tnmom1c(0, x, +Inf, 0, 1) ≈ TN.tnmean(x, +Inf)
    @test TN.tnmom1c(0, -Inf, x, 0, 1) ≈ TN.tnmean(-Inf, x)
end

@test TN.tnmom1c(10, 0, 10) ≈ -9.20211543919713464412024961258
@test TN.tnmom1c(1, 3, 4) ≈ 2.26045428559002115390313889722
@test TN.tnmom1c(-1, 3, 4) ≈ 4.26045428559002115390313889722
@test TN.tnmom1c(2, 3, 4) ≈ 1.26045428559002115390313889722
@test TN.tnmom1c(-2, 3, 4) ≈ 5.26045428559002115390313889722
@test TN.tnmom1c(0, -1e200, 1e200) ≈ 0
@test TN.tnmom1c(1, -1e200, 1e200) ≈ -1
@test TN.tnmom1c(1e200, -1e200, 1e200) ≈ -1e200
@test TN.tnmom1c(0, -1e6, -999000) ≈ TN.tnmean(-1e6, -999000)
@test TN.tnmom1c(0, +1e6, +Inf) ≈ TN.tnmean(+1e6, +Inf)
@test TN.tnmom1c(0, -Inf, -1e6) ≈ TN.tnmean(-Inf, -1e6)

@test TN.tnmom1c(999000, 999000, 1e6) ≈ 1.00100100099899498898098100910e-6
@test TN.tnmom1c(999500, 999000, 1e6) ≈ -499.999998998998999001005011019
@test TN.tnmom1c(1e6, 999000, 1e6) ≈ -999.999998998998999001005011019
