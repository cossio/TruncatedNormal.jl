using Statistics, SpecialFunctions
export tnvar

"""
    tnvar(a, b)

Variance of the standard normal distribution truncated to the interval [a,b].
"""
function tnvar(a, b)
    if !(a ≤ b)
        return oftype(middle(a, b), NaN)
    elseif a == b
        return zero(middle(a, b))
    elseif abs(a) > abs(b)
        return tnvar(-b, -a)
    end

    @assert a < b && abs(a) ≤ abs(b)
    @assert a ≤ 0 ≤ b || 0 ≤ a < b
    
    return tnmom2(a, b) - tnmom1(a, b)^2
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
