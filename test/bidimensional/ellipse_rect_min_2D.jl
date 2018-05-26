using Base.Test

import TruncatedNormal: ellipse_rect_min_2D

@testset "ellipse_rect_min_2D" begin
    srand(0)
    for repetition = 1 : 1000
        ρ = 2rand() - 1
        a1, b1 = minmax((10rand(2) - 5)...)
        a2, b2 = minmax((10rand(2) - 5)...)
        if -1 < ρ < 1
            x1, x2 = ellipse_rect_min_2D(ρ, (a1, a2), (b1, b2))
            y1 = a1 + (b1 - a1) * rand()
            y2 = a2 + (b2 - a2) * rand()
            #@show ρ, (a1,a2), (b1,b2)
            #@show (x1,x2), (y1,y2)
            @test x1^2 + x2^2 - 2ρ * x1 * x2 ≤ y1^2 + y2^2 - 2ρ * y1 * y2
        end
    end
end
