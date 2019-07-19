using SafeTestsets

@safetestset "tnmom1" begin include("tnmom1.jl") end
@safetestset "tnmom2" begin include("tnmom2.jl") end
@safetestset "tnvar" begin include("tnvar.jl") end
