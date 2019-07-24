using Statistics, SpecialFunctions
export tnmom1, tnmean

"""
    tnmom1(a,b)

Mean of the truncated standard normal distribution.
"""
function tnmom1(a::Real, b::Real)
    if !(a ≤ b)
        return oftype(middle(a, b), NaN)
    elseif a == b
        return middle(a, b)
    elseif abs(a) > abs(b)
        return -tnmom1(-b, -a)
    elseif isinf(a) && isinf(b)
        return zero(middle(a, b))
    end

    @assert a < b && abs(a) ≤ abs(b)
    @assert a ≤ 0 ≤ b || 0 ≤ a < b

    Δm1 = expm1((a - b)middle(a, b))
    Δ = one(Δm1) + Δm1

    if 0 ≤ a < b && b ≥ 1
        return √(2/π) * Δm1 / (Δ * erfcx(b / √2) - erfcx(a / √2))
    else # a ≤ 0 ≤ b
        return √(2/π) * Δm1 * exp(-a^2 / 2) / (erf(a / √2) - erf(b / √2))
    end
end

"""
    tnmom1(a, b, μ, σ)

Mean of the truncated normal distribution, where μ, σ are the mean and standard
deviation of the untruncated distribution.
"""
function tnmom1(a, b, μ, σ)
    α = (a - μ) / σ
    β = (b - μ) / σ
    return μ + tnmom1(α, β) * σ
end

"""
    tnmean(a, b)

Mean of the truncated standard normal distribution.

    tnmean(a, b, μ, σ)

Mean of the truncated normal distribution, where μ, σ are the mean and standard
deviation of the untruncated distribution.
"""
tnmean = tnmom1
