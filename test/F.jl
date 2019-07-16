using Test, SpecialFunctions
using TruncatedNormal: _F1, _F2, _F3

f1(x, y) = (exp(-x^2) - exp(-y^2)) / (erf(x) - erf(y))
f2(x, y) = (exp(-x^2)x - exp(-y^2)y) / (erf(x) - erf(y))
f3(x, y) = (exp(-x^2) + exp(-y^2)) / (erf(x) - erf(y))

@testset "F1, F2, F3" begin
    @testset "x ≈ y" begin
        @test _F1(1, 1 + 1e-8) ≈ -1.77245385976778525228484995704
        @test _F2(1, 1 + 1e-8) ≈ -0.886226943177296522704243423834
        @test _F3(1, 1 + 1e-8) ≈ -1.77245385090551605683906558925e8
        @test _F1(0, 1e-10) ≈ -8.86226925452758013647606696795e-11
        @test _F2(0, 1e-10) ≈ 0.886226925452758013643175562168
        @test _F3(0, 1e-10) ≈ -1.77245385090551602729521339359e10
    end

    @testset "Properties" begin
        for x = -10 : 1e-2 : 10
            @test _F1(x, x) ≈ -√π * x
            @test _F2(x, x) ≈ √π * (1 - 2x^2) / 2
            @test iszero(_F1(x, -x))
            if x ≠ 0
                @test _F2(x, -x) ≈ x * exp(-x^2) / erf(x)
                @test _F3(x, -x) ≈ exp(-x^2) / erf(x)
            end
        end
        for x = -10 : 0.1 : 10, y = -10 : 0.1 : 10
            @test _F1(x, y) == _F1(y, x)
            @test _F2(x, y) == _F2(y, x)
            @test _F1(-x, -y) == -_F1(x, y)
            @test _F2(-x, -y) == _F2(x, y)
            if x ≠ y
                @test _F3(x, y) == -_F3(y, x)
                @test _F3(-x, -y) == -_F3(x, y)
            end
        end
        for x = -1e-8 : 1e-9 : 1e-8
            if x ≠ 0
                @test _F2(x, -x) ≈ x * exp(-x^2) / erf(x)
                @test _F3(x, -x) ≈ exp(-x^2) / erf(x)
            end
        end
    end

    @testset "x ≤ 0 ≤ y || y ≤ 0 ≤ x" begin
        @test _F1(-1, 1) == 0
        @test _F2(-1, 1) ≈ 0.436548113220292413451721744841
        @test _F3(-1, 1) ≈ -0.436548113220292413451721744841
        @test _F1(-2, 1) ≈ 0.190184666491092019086029997367
        @test _F2(-2, 1) ≈ 0.220079240679365999668879344802
        @test _F1(-100, 1) ≈ 0.199641440747717373738835579396
        @test _F2(-100, 1) ≈ 0.199641440747717373738835579396
        @test _F3(-100, 1) ≈ -0.199641440747717373738835579396
        @test _F1(-1, 100) ≈ -0.199641440747717373738835579396
        @test _F2(-1, 100) ≈ 0.199641440747717373738835579396
        @test _F3(-1, 100) ≈ -0.199641440747717373738835579396

        @test _F1(-Inf, 0) ≈ 1
        @test _F1(-Inf, 1) ≈ 0.199641440747717373738835579396
        @test _F1(-Inf, 4) ≈ 5.62675877933755134754466836019e-8
        @test _F1(-Inf, 10) ≈ 1.86003798801041798147984790193e-44
    end

    @testset "x, y < 0" begin
        @test _F1(-2, -1) ≈ 2.29039726549175154756498826419
        @test _F2(-2, -1) ≈ -2.17039030552464315391802895425
        @test _F1(-101, -100) ≈ 177.254246473800679651649984682
        @test _F2(-101, -100) ≈ -17725.424647380067965164998468
        @test _F3(-101, -100) ≈ -177.254246473800679651649984683
        @test _F1(-110, -150) ≈ 194.977979542322091092377246367
        @test _F2(-110, -150) ≈ -21447.577749655430020161497100
        @test _F3(-110, -150) ≈ 194.977979542322091092377246367
    end

    @testset "x, y > 0" begin
        @test _F1(1, 2) ≈ -2.29039726549175154756498826419
        @test _F2(1, 2) ≈ -2.17039030552464315391802895425
        @test _F1(100, 101) ≈ -177.254246473800679651649984683
        @test _F2(100, 101) ≈ -17725.4246473800679651649984683
        @test _F1(110, 150) ≈ -194.977979542322091092377246367
        @test _F2(110, 150) ≈ -21447.5777496554300201614971004
        @test _F1(100, 115) ≈ -177.254246473800679651649984683
        @test _F2(100, 115) ≈ -17725.4246473800679651649984683
    end

    @testset "Exhaustive" begin
        for x = (-big"5e-8", big"5e-8")
            @test _F1(Float64(x), Float64(-x)) ≈ f1(x, -x)
            @test _F2(Float64(x), Float64(-x)) ≈ f2(x, -x)
            @test _F3(Float64(x), Float64(-x)) ≈ f3(x, -x)
        end
        for x = -1 : 0.01 : 1, y = -1 : 0.01 : 1
            x == y && continue
            @test _F1(x, y) ≈ f1(x, y)
            @test _F2(x, y) ≈ f2(x, y)
            @test _F3(x, y) ≈ f3(x, y)
        end
    end
end


@testset "F1 & F2 (infinite arguments)" begin
    @test _F1(Inf, Inf) == -Inf
    @test _F1(-Inf, -Inf) == Inf
    @test _F1(-Inf, Inf) == _F1(Inf, -Inf) == 0

    @test _F2(-Inf, Inf) == _F2(Inf, -Inf) == 0
    @test _F2(Inf, Inf) == _F2(-Inf, -Inf) == -Inf

    for v in (0, 1, -1)
        @test _F1(-Inf, v) ≈ _F1(v, -Inf) ≈ _F1(-1e20, v) ≈ _F1(v, -1e20)
        @test _F2(-Inf, v) ≈ _F2(v, -Inf) ≈ _F2(-1e20, v) ≈ _F2(v, -1e20)
        @test _F1(Inf, v) ≈ _F1(v, Inf) ≈ _F1(1e20, v) ≈ _F1(v, 1e20)
        @test _F2(Inf, v) ≈ _F2(v, Inf) ≈ _F2(1e20, v) ≈ _F2(v, 1e20)
    end

    for v = -10 : 0.1 : 10
        @test _F1(-Inf, v) ≈ _F1(v, -Inf) ≈ _F1(-1e20, v) ≈ _F1(v, -1e20)
        @test _F2(-Inf, v) ≈ _F2(v, -Inf) ≈ _F2(-1e20, v) ≈ _F2(v, -1e20)
        @test _F1(Inf, v) ≈ _F1(v, Inf) ≈ _F1(1e20, v) ≈ _F1(v, 1e20)
        @test _F2(Inf, v) ≈ _F2(v, Inf) ≈ _F2(1e20, v) ≈ _F2(v, 1e20)
    end
end