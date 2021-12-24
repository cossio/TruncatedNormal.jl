using Test, SpecialFunctions
import TruncatedNormal as TN

@test TN.tnmean(-Inf, Inf) == 0
@test TN.tnmean(0, +Inf) ≈ +√(2/π)
@test TN.tnmean(-Inf, 0) ≈ -√(2/π)
@test TN.tnmean(-Inf, -Inf) == -Inf
@test TN.tnmean(+Inf, +Inf) == +Inf

for x = -10:10
    @test TN.tnmean(x, +Inf) ≈ +√(2/π) / erfcx(+x/√2)
    @test TN.tnmean(-Inf, x) ≈ -√(2/π) / erfcx(-x/√2)
end

for x = exp.(-100:100)
    @test TN.tnmean(-x, x) ≈ 0
    @test TN.tnmean(0, x) ≈ -√(2/π) * expm1(-x^2 / 2) / erf(x / √2)
end

@test TN.tnmean(1e-44, 1e-43) ≈ 5.4999999999999999999999999999999999999999e-44

@test TN.tnmean(-Inf, Inf, 0, 1) == TN.tnmean(-Inf, Inf)
for x = -100:100
    @test TN.tnmean(x, +Inf, 0, 1) ≈ TN.tnmean(x, +Inf)
    @test TN.tnmean(-Inf, x, 0, 1) ≈ TN.tnmean(-Inf, x)
end

@test TN.tnmean(-1e4, 1e4, 0, 1) == TN.tnmean(-1e4, 1e4)

for x = -100:100, y = x + 1 : 100
    @test TN.tnmean(x, y, 0, 1) == TN.tnmean(x, y)
end

@test TN.tnmean(100, 115) ≈ 100.00999800099926070518490239457545847490332879043
@test TN.tnmean(-1e6, -999000) ≈ -999000.00000100100100099899498898098
@test TN.tnmean(+1e6, +Inf) ≈ +1.00000000000099999999999800000e6
@test TN.tnmean(-Inf, -1e6) ≈ -1.00000000000099999999999800000e6

@test TN.tnmean(-1e200, 1e200) ≈ 0
@test TN.tnmean(0, +1e200) ≈ +0.797884560802865355879892119869
@test TN.tnmean(-1e200, 0) ≈ -0.797884560802865355879892119869

@test TN.tnmean(50, 70, -2, 3) ≈ 50.171943499898757645751683644632860837133138152489
@test TN.tnmean(-100.0, 0.0, 0.0, 2.0986317998643735) ≈ -1.6744659119217125058885983754999713622460154892645
@test TN.tnmean(0.0, 0.9, 0.0, 0.07132755843183151) ≈ 0.056911157632522598806524588414964004271754161737065
@test TN.tnmean(-100.0, 100.0, 0.0, 17.185261847875548) ≈ 0
@test TN.tnmean(-100.0, 0.5, 0.0, 0.47383322897860064) ≈ -0.1267981330521791493635176736743283314399
@test TN.tnmean(-100.0, 100.0, 0.0, 17.185261847875548) ≈ 0

for x = 1e1 .^ (-400:100:400)
    @test TN.tnmean(-x, x) == TN.tnmean(-x, x) || isnan(TN.tnmean(-x, x)) && isnan(TN.tnmean(-x, x))
end

for a = exp.(-10:10), b = exp.(-10:10)
    a ≤ b || continue
    @test a ≤ TN.tnmean(a, b) ≤ b
end

# https://github.com/JuliaStats/Distributions.jl/issues/827
@test TN.tnmean(0, 1000, 1000000, 1) ≈ 999.99999899899899900100501101901899090472046236710608108591983
