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
    
    if a ≤ 0 ≤ b || -1 ≤ a ≤ b ≤ 1
    elseif 0 ≤ a < b
        #m1 = 
        m2 = 2/√π * _F2(a/√2, b/√2)
        m1 = √(2/π) * _F1(a/√2, b/√2)
        return 1 - (middle(a,b) - (b-a)/2 * _G3(a/√2, b/√2) + m1) * m1
    # elseif 0 ≤ a ≤ b
    #     c2 = 1 + a^2 - 2/√π * _F4(a/√2, b/√2)
    #     return c2 - (m1 - a)^2
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
