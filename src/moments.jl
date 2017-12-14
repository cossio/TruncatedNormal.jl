export tnmean, tnvar

"""
    tnmean(a,b)
Mean of the truncated standard normal distribution

    tnmean(a, b, μ, σ)
Mean of the truncated normal distribution, where μ,σ
are the mean and standard deviation of the untruncated
distribution.
"""
function tnmean end

function tnmean(a, b)
    -Inf < a ≤ b < Inf || throw(ArgumentError("tnmean called with invalid parameters a=$a, b=$b"))
    √(2/π) * _F1(a/√2, b/√2)
end

function tnmean(a, b, μ, σ)
    -Inf < a ≤ b < Inf || throw(ArgumentError("tnmean called with invalid parameters a=$a, b=$b"))
    -Inf < μ < Inf && 0 < σ < Inf || throw(DomainError())
    α = (a - μ)/σ; β = (b - μ)/σ
    μ + tnmean(α, β) * σ
end


"""
    tnvar(a,b)
Variance of the truncated standard normal distribution

    tnvar(a, b, μ, σ)
Variance of the truncated normal distribution, where μ,σ
are the mean and standard deviation of the untruncated
distribution.
"""
function tnvar end

function tnvar(a, b)
    -Inf < a ≤ b < Inf || throw(ArgumentError("tnvar called with invalid parameters a=$a, b=$b"))
    1 + 2/√π * _F2(a/√2, b/√2) - 2/π * _F1(a/√2, b/√2)^2
end

function tnvar(a, b, μ, σ)
    -Inf < a ≤ b < Inf || throw(ArgumentError("tnvar called with invalid parameters a=$a, b=$b"))
    -Inf < μ < Inf && 0 < σ < Inf || throw(DomainError())
    α = (a - μ)/σ; β = (b - μ)/σ
    tnvar(α, β) * σ^2
end
