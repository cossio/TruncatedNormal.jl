using Base.Test

import TruncatedNormal
TN = TruncatedNormal

@testset "moments" begin
    @test TN.tnmean(100, 115) ≈ 100.00999800099926070518490239457545847490332879043
    
end
