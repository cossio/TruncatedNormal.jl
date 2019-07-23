using Test, TruncatedNormal, SpecialFunctions

@test tnmom1(-Inf, Inf) == 0
@test tnmom1(0, +Inf) ≈ +√(2/π)
@test tnmom1(-Inf, 0) ≈ -√(2/π)
@test tnmom1(-Inf, -Inf) == -Inf
@test tnmom1(+Inf, +Inf) == +Inf

for x = -10:10
    @test tnmom1(x, +Inf) ≈ +√(2/π) / erfcx(+x/√2)
    @test tnmom1(-Inf, x) ≈ -√(2/π) / erfcx(-x/√2)
end

for x = exp.(-100:100)
    @test tnmom1(-x, x) ≈ 0
end

@test tnmom1(-Inf, Inf, 0, 1) == tnmom1(-Inf, Inf)

for x = -100:100, y = x + 1 : 100
    @test tnmom1(x, y, 0, 1) == tnmom1(x, y)
end

for x = -100:100
    @test tnmom1(x, +Inf, 0, 1) ≈ tnmom1(x, +Inf)
    @test tnmom1(-Inf, x, 0, 1) ≈ tnmom1(-Inf, x)
end

@test tnmom1(-1e4, 1e4, 0, 1) == tnmom1(-1e4, 1e4)

@test tnmom1(100, 115) ≈ 100.00999800099926070518490239457545847490332879043
@test tnmom1(-1e6, -999000) ≈ -999000.00000100100100099899498898098
@test tnmom1(+1e6, +Inf) ≈ +1.00000000000099999999999800000e6
@test tnmom1(-Inf, -1e6) ≈ -1.00000000000099999999999800000e6

@test tnmom1(-1e200, 1e200) ≈ 0
@test tnmom1(0, +1e200) ≈ +0.797884560802865355879892119869
@test tnmom1(-1e200, 0) ≈ -0.797884560802865355879892119869

@test tnmom1(50, 70, -2, 3) ≈ 50.171943499898757645751683644632860837133138152489
@test tnmom1(-100.0, 0.0, 0.0, 2.0986317998643735) ≈ -1.6744659119217125058885983754999713622460154892645
@test tnmom1(0.0, 0.9, 0.0, 0.07132755843183151) ≈ 0.056911157632522598806524588414964004271754161737065
@test tnmom1(-100.0, 100.0, 0.0, 17.185261847875548) ≈ 0
@test tnmom1(-100.0, 0.5, 0.0, 0.47383322897860064) ≈ -0.1267981330521791493635176736743283314399
@test tnmom1(-100.0, 100.0, 0.0, 17.185261847875548) ≈ 0

for x = 1e1 .^ (-400:100:400)
    @test tnmean(-x, x) == tnmom1(-x, x) || isnan(tnmean(-x, x)) && isnan(tnmom1(-x, x))
end
