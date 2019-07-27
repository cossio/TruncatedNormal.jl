using Statistics, SpecialFunctions
export tnvar

"""
    tnvar(a, b)

Variance of the standard normal distribution truncated to the interval [a,b].
"""
function tnvar(a::Real, b::Real)
    if a == b
        return zero(middle(a, b))
    elseif a < b
        m1 = tnmom1(a, b)
        m2 = √tnmom2(a, b)
        return (m2 - m1) * (m2 + m1)

        m1 = tnmom1(a, b)
        @assert a ≤ m1 ≤ b
        return tnmom2c(m1, a, b)
    else
        return oftype(middle(a, b), NaN)
    end
end

"""
    tnvar(a, b, μ, σ)

Variance of the truncated normal distribution, where μ, σ are the mean and
standard deviation of the untruncated distribution.
"""
function tnvar(a, b, μ, σ)
    α = (a - μ) / σ
    β = (b - μ) / σ
    return tnvar(α, β) * σ^2
end
