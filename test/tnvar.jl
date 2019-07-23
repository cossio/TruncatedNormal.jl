using Test, TruncatedNormal, SpecialFunctions, Statistics

@test isnan(tnvar(1, 0))
@test isnan(tnvar(+Inf, 0))
@test isnan(tnvar(0, -Inf))
@test isnan(tnvar(NaN, NaN))
@test isnan(tnvar(0, NaN))
@test isnan(tnvar(NaN, 0))

@test tnvar(-Inf, +Inf) == 1
@test tnvar(+Inf, +Inf) == 0
@test tnvar(-Inf, -Inf) == 0
@test tnvar(0, +Inf) ≈ 1 - 2/π
@test tnvar(-Inf, 0) ≈ 1 - 2/π

for x = -10:10
    @test tnvar(x, +Inf) ≈ 1 + √(2/π) * x / erfcx(+x / √2) - 2/π / erfcx(+x / √2)^2
    @test tnvar(-Inf, x) ≈ 1 - √(2/π) * x / erfcx(-x / √2) - 2/π / erfcx(-x / √2)^2
end

for x = -10:10
    @test tnvar(x, +Inf, 0, 1) ≈ tnvar(x, +Inf)
    @test tnvar(-Inf, x, 0, 1) ≈ tnvar(-Inf, x)
end

@test tnvar(50, 70) ≈ 0.0003990431868038995479099272265360593305365
@test tnvar(1e6, Inf) ≈ 9.99999999994000000000050000000e-13

# https://github.com/JuliaStats/Distributions.jl/issues/827
@test_broken tnvar(999000, 1e6) ≥ 0
@test_broken tnvar(999000, 1e6) ≈ 1.00200300399898194688784897455e-12

@test tnvar(-Inf, Inf, 0, 1) == 1

@test tnvar(50, 70, -2, 3) ≈ 0.029373438107168350377591231295634273607812172191712
@test tnvar(-100.0, 0.0, 0.0, 2.0986317998643735) ≈ 1.6004193412141677189841357987638847137391508803335
@test tnvar(0.0, 0.9, 0.0, 0.07132755843183151) ≈ 0.0018487407287725028827020557707636415445504260892486
@test tnvar(-100.0, 100.0, 0.0, 17.185261847875548) ≈ 295.333163899557735486302841237124507431445
@test tnvar(-100.0, 0.5, 0.0, 0.47383322897860064) ≈ 0.145041095812679283837328561547251019229612
@test tnvar(-100.0, 100.0, 0.0, 17.185261847875548) ≈ 295.333163899557735486302841237124507431445

@test tnvar(-10000, 10000, 0, 1) ≈ 1
