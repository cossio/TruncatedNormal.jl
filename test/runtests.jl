using SafeTestsets

@safetestset "util" begin include("util.jl") end
@safetestset "tnmean" begin include("tnmean.jl") end
@safetestset "tnmom2" begin include("tnmom2.jl") end
@safetestset "tnmom1c" begin include("tnmom1c.jl") end
@safetestset "tnmom2c" begin include("tnmom2c.jl") end
@safetestset "tnvar" begin include("tnvar.jl") end
@safetestset "truncnorm" begin include("truncnorm.jl") end
@safetestset "rejection" begin include("rejection.jl") end
