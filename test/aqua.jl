import TruncatedNormal
import Aqua

using Test: @testset

@testset "Aqua" begin
    Aqua.test_all(TruncatedNormal)
end
