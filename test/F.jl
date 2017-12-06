using Base.Test

import TruncatedNormal
TN = TruncatedNormal

@testset "F1 & F2" begin
    @testset "|x-y| ≤ thresh" begin
        @test TN.F1(1, 1 + 1e-8; thresh=1e-7) ≈ 1.7724538597677852522848499570355616296194525129226
        @test TN.F2(1, 1 + 1e-8; thresh=1e-7) ≈ 0.88622694317729652270424342383429653215361427389442
    end

    @testset "x ≤ 0 ≤ y || y ≤ 0 ≤ x" begin
        @test TN.F1(-1,1) == 0
        @test TN.F1(-2,1) ≈ -0.19018466649109201908602999736689089421341283247396
        @test TN.F2(-1,1) ≈ -0.43654811322029241345172174484122488428321404645080
        @test TN.F1(-100,1) ≈ -0.19964144074771737373883557939628614775906698169272
        @test TN.F2(-100,1) ≈ -0.19964144074771737373883557939628614775906698169272
        @test TN.F1(-1,100) ≈ 0.19964144074771737373883557939628614775906698169272
        @test TN.F2(-1,100) ≈ -0.19964144074771737373883557939628614775906698169272
    end

    @testset "x,y < 0" begin
        @test TN.F1(-2,-1) ≈ -2.2903972654917515475649882641880377756080219194918
        @test TN.F2(-2,-1) ≈ 2.1703903055246431539180289542483966628176105020029
        @test TN.F1(-101,-100) ≈ -177.25424647380067965164998468281027353208843866970
        @test TN.F2(-101,-100) ≈ 17725.424647380067965164998468281027353208843866970
        @test TN.F1(-110,-150) ≈ -194.97797954232209109237724636721466012249482102215
        @test TN.F2(-110,-150) ≈ 21447.577749655430020161497100393612613474430312437
    end

    @testset "x,y > 0" begin
        @test TN.F1(1,2) ≈ 2.2903972654917515475649882641880377756080219194918
        @test TN.F2(1,2) ≈ 2.1703903055246431539180289542483966628176105020029
        @test TN.F1(100,101) ≈ 177.25424647380067965164998468281027353208843866970
        @test TN.F2(100,101) ≈ 17725.424647380067965164998468281027353208843866970
        @test TN.F1(110,150) ≈ 194.97797954232209109237724636721466012249482102215
        @test TN.F2(110,150) ≈ 21447.577749655430020161497100393612613474430312437

        @test TN.F1(100,115) ≈ 177.25424647380067965164998468281027353208843866970
        @test TN.F2(100,115) ≈ 17725.424647380067965164998468281027353208843866970
    end

    @testset "F(x,y) == F(y,x)" begin
        for r = 1:10
            x,y = rand(2) - 0.5
            tests = [(@test TN.F1(x,y) == TN.F1(y,x)),
                     (@test TN.F2(x,y) == TN.F2(y,x))]
            any(typeof.(tests) .≠ Base.Test.Pass) && @show x,y
        end
    end
end