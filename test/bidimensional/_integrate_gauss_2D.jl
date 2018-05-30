using Base.Test

import TruncatedNormal: _integrate_gauss_2D

@testset "_integrate_gauss_2D" begin
    v, C = _integrate_gauss_2D(x -> 1., 0.5, (-0.5,0.2), (-0.2,1.)); @test v * exp(C) ≈ 0.02752462487030073;
    v, C = _integrate_gauss_2D(x -> 1., 0.5, (0.5, 0.8), (0.8, 1.)); @test v * exp(C) ≈ 0.007118405893210067;
    v, C = _integrate_gauss_2D(x -> x[1]*x[2], 0.5, (0.5, 0.8), (0.8, 1.)); @test v * exp(C) ≈ 0.004139760367882591;
    v, C = _integrate_gauss_2D(x -> x[1], 0.5, (0.5, 0.8), (0.8, 1.)); @test v * exp(C) ≈ 0.004612694953630602;
end
