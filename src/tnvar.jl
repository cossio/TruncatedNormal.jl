using Statistics, SpecialFunctions
export tnvar

"""
    tnvar(a, b)

Variance of the standard normal distribution truncated to the interval [a,b].
"""
function tnvar(a::Real, b::Real)
    if !(a ≤ b)
        return oftype(middle(a, b), NaN)
    elseif a == b
        return zero(middle(a, b))
    else
        m1 = tnmom1(a, b)
        m2 = √tnmom2(a, b)
        return (m2 - m1) * (m2 + m1)

        m = tnmom1(a, b)
        return tnmom2c(m, a, b)
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
