#=
Functions to compute the sample the truncated bivariate Gaussian.
Uses an R library.
=#

module R

import RCall

RealMatrix = AbstractMatrix{<:Real}
RealVector = AbstractVector{<:Real}

"""
Samples N points from a bivariate Gaussian X ~ N(0,Sig) 
conditional on lb ≤ X ≤ ub. Accepts infinite lb, ub.
"""
function sample_truncated_bivariate_gaussian(lb::RealVector, ub::RealVector, Sig::RealMatrix, N::Integer)
    #=
    This makes a call to the Z. I. Botev's R package TruncatedNormal.
    =#
    @assert size(Sig) == (2,2) && length(lb) == length(ub) == 2
    @assert N > 0
    @assert all(lb .< ub)
    @assert issymmetric(Sig)
    
    RCall.R"library(TruncatedNormal)"

    RCall.@rput lb
    RCall.@rput ub
    RCall.@rput Sig
    RCall.@rput N

    RCall.R"X <- TruncatedNormal::mvrandn(l = lb, u = ub, Sig = Sig, n = N)"
    RCall.@rget X

    return X
end


"""
Samples N points from a bivariate Gaussian X ~ N(μ,Sig) 
conditional on lb ≤ X ≤ ub. Accepts infinite lb, ub.
"""
function sample_truncated_bivariate_gaussian(lb::RealVector, ub::RealVector, μ::RealVector, Sig::RealMatrix, N::Integer)
    @assert size(Sig) == (2,2) && length(μ) == length(lb) == length(ub) == 2
    @assert N > 0
    @assert all(lb .< ub)
    @assert all(-Inf .< μ .< Inf)
    @assert issymmetric(Sig)

    X = sample_truncated_bivariate_gaussian(lb .- μ, ub .- μ, Sig, N)
    return X .+ μ
end


"""
Compute the mean and covariance matrix of a Gaussian
bivariate X ~ N(μ,Σ), conditional on lb < X < ub.
"""
function truncated_bivariate_gaussian_moments(lb::RealVector, ub::RealVector, μ::RealVector, Σ::RealMatrix, N::Integer = 10000)
    X = sample_truncated_bivariate_gaussian(lb, ub, μ, Σ, N)
    tμ = mean(X, 2)
    tΣ = [mean(X[i,:], X[j,:]) for i = 1:2, j = 1:2]
    return tμ, tΣ
end


end