using Base.Test

import TruncatedNormal: gauss2Dmoment1, gauss2Dmoment2, gauss2Dmoment11, gauss2Dmoment22, gauss2Dmoment12

@testset "gauss2Dmoments" begin
    a = (-0.5, 0.2); b = (-0.2, 1.); 
    
    μ = (0,0); Σ = [1 0.5; 0.5 1]
    @test gauss2Dmoment1(μ, Σ, a, b) ≈ -0.3437958060623487
    @test gauss2Dmoment2(μ, Σ, a, b) ≈ 0.5472225561505393
    @test gauss2Dmoment11(μ, Σ, a, b) ≈ 0.12564309156273065
    @test gauss2Dmoment22(μ, Σ, a, b) ≈ 0.34966461022240974
    @test gauss2Dmoment12(μ, Σ, a, b) ≈ -0.18788358600082938

    μ = (0, 0); Σ = [1.2 -0.1; -0.1 2.3]
    @test gauss2Dmoment1(μ, Σ, a, b) ≈ -0.3479699198765843
    @test gauss2Dmoment2(μ, Σ, a, b) ≈ 0.5868430723295551
    @test gauss2Dmoment11(μ, Σ, a, b) ≈ 0.12856179814507035
    @test gauss2Dmoment22(μ, Σ, a, b) ≈ 0.3971203030711006
    @test gauss2Dmoment12(μ, Σ, a, b) ≈ -0.20421807848632167

    μ = (0.2, -0.1); Σ = [1 -0.1; -0.1 1]
    @test gauss2Dmoment1(μ, Σ, a, b) ≈ -0.34635036544024095
    @test gauss2Dmoment2(μ, Σ, a, b) ≈ 0.5661235115386849
    @test gauss2Dmoment11(μ, Σ, a, b) ≈ 0.12742791361156347
    @test gauss2Dmoment22(μ, Σ, a, b) ≈ 0.37201427804705994
    @test gauss2Dmoment12(μ, Σ, a, b) ≈ -0.19611595384101185

    μ = (0.2, -0.1); Σ = [1.2 -0.1; -0.1 2.3]
    @test gauss2Dmoment1(μ, Σ, a, b) ≈ -0.34674589834822167
    @test gauss2Dmoment2(μ, Σ, a, b) ≈ 0.5849236447178073
    @test gauss2Dmoment11(μ, Σ, a, b) ≈ 0.12770757916520126
    @test gauss2Dmoment22(μ, Σ, a, b) ≈ 0.39483894539106945
    @test gauss2Dmoment12(μ, Σ, a, b) ≈ -0.20283420009833092
end


@testset "gauss2Dmoments, far into the tails" begin
    a = (100, 200); b = (120, 230); 
    μ = (0,0); Σ = [1 0.9; 0.9 1];
    # the following values were obtained using matlab/sample_bivariate_trunc.m
    @test gauss2Dmoment1(μ, Σ, a, b) ≈ 1.199968324383764e2  rtol=1e-7
    @test gauss2Dmoment2(μ, Σ, a, b) ≈ 2.000020659528939e2  rtol=1e-7
    @test gauss2Dmoment11(μ, Σ, a, b) ≈ 1.439923980527374e4 rtol=1e-7
    @test gauss2Dmoment22(μ, Σ, a, b) ≈ 4.000082638968676e4 rtol=1e-7
    @test gauss2Dmoment12(μ, Σ, a, b) ≈ 2.399961439547430e4 rtol=1e-7
end
