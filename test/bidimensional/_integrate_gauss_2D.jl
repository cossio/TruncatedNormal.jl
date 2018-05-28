using Base.Test

import TruncatedNormal: _integrate_gauss_2D

@testset "_integrate_gauss_2D" begin
    @test _integrate_gauss_2D(x -> 1., 0.5, (-0.5,0.2), (-0.2,1.), 0) ≈ 0.02752462487030073
    @test _integrate_gauss_2D(x -> 1., 0.5, (0.5, 0.8), (0.8, 1.), 0) ≈ 0.007118405893210067
    @test _integrate_gauss_2D(x -> x[1]*x[2], 0.5, (0.5, 0.8), (0.8, 1.), 0) ≈ 0.004139760367882591
    @test _integrate_gauss_2D(x -> x[1], 0.5, (0.5, 0.8), (0.8, 1.), 0) ≈ 0.004612694953630602
end
