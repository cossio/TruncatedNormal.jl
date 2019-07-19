using Statistics, SpecialFunctions
export tnmom1, tnmom2, tnvar, tnmean

"""
    tnmom2(a, b)

Second moment of the truncated standard normal distribution.
"""
function tnmom2(a, b)
    if !(a ≤ b)
        return oftype(middle(a, b), NaN)
    elseif a == b
        return oftype(middle(a, b), a * b)
    elseif abs(a) > abs(b)
        return tnmom2(-b, -a)
    elseif isinf(a) && isinf(b)
        return one(middle(a, b))
    elseif isfinite(a) && isinf(b)
        return 1 + √(2 / π) * a / erfcx(a / √2)
    end

    @assert a < b < Inf && abs(a) ≤ abs(b)
    @assert a ≤ 0 ≤ b || 0 ≤ a ≤ b

    Δ = exp((a - b)middle(a, b))
    if a ≤ 0 ≤ b
        return 1 - √(2 / π) * (Δ * b - a)exp(-a^2 / 2) / (erf(b / √2) - erf(a / √2))
    elseif 0 ≤ a ≤ b
        return 1 - √(2 / π) * (Δ * b - a) / (erfcx(a / √2) - Δ * erfcx(√2))
    end
end

"""
    tnmom2(a, b, μ, σ)

Second moment of the truncated normal distribution, where μ, σ
are the mean and standard deviation of the untruncated distribution.
"""
tnmom2(a, b, μ, σ) = tnmom1(a, b, μ, σ)^2 + tnvar(a, b, μ, σ)
