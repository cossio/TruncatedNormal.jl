using Base.Test

import TruncatedNormal: _ellipse_rect_min_2D, _ellipse_rect_max_2D, _ellipse_quadratic

@testset "ellipse_rect_min_2D" begin
    srand(2)
    for repetition = 1 : 1000
        ρ = 2rand() - 1
        a1, b1 = minmax((10rand(2) - 5)...)
        a2, b2 = minmax((10rand(2) - 5)...)
        if -1 < ρ < 1
            x1, x2 = _ellipse_rect_min_2D(ρ, (a1, a2), (b1, b2))
            y1 = a1 + (b1 - a1) * rand()
            y2 = a2 + (b2 - a2) * rand()
            @test _ellipse_quadratic(ρ, (x1, x2)) ≤ _ellipse_quadratic(ρ, (y1, y2))
            #@test x1^2 + x2^2 - 2ρ * x1 * x2 ≤ y1^2 + y2^2 - 2ρ * y1 * y2
        end
    end
end


@testset "ellipse_rect_max_2D" begin
    srand(2)
    for repetition = 1 : 1000
        ρ = 2rand() - 1
        a1, b1 = minmax((10rand(2) - 5)...)
        a2, b2 = minmax((10rand(2) - 5)...)
        if -1 < ρ < 1
            x1, x2 = _ellipse_rect_max_2D(ρ, (a1, a2), (b1, b2))
            y1 = a1 + (b1 - a1) * rand()
            y2 = a2 + (b2 - a2) * rand()
            @test _ellipse_quadratic(ρ, (x1, x2)) ≥ _ellipse_quadratic(ρ, (y1, y2))
            #@test x1^2 + x2^2 - 2ρ * x1 * x2 ≥ y1^2 + y2^2 - 2ρ * y1 * y2
        end
    end
end
