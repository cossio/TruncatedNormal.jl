using Statistics, SpecialFunctions
export tnmom1, tnmean, tnmom1i

"""
    tnmom1(a, b)

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
    @assert a ≤ 0 ≤ b || 0 < a < b

    Δ = (b - a) * middle(a, b)
    #Δ = one(Δm1) + Δm1

    if a ≤ 0 ≤ b
        m = √(2/π) * expm1(-Δ) * exp(-a^2 / 2) / erf(b/√2, a/√2)
    elseif 0 < a < b
        z = exp(-Δ) * erfcx(b/√2) - erfcx(a/√2)
        iszero(z) && return middle(a, b)
        m = √(2/π) * expm1(-Δ) / z
    end
    return clamp(m, a, b)
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

"""
    tnmom1i(a, b)

Mean of the normal distribution with variance -1 and mean 0, truncated to [a,b].
"""
function tnmom1i(a::Real, b::Real)
    if !(-Inf < a ≤ b < Inf)
        return oftype(middle(a, b), NaN)
    elseif a == b
        return middle(a, b)
    elseif abs(a) > abs(b)
        return -tnmom1i(-b, -a)
    end

    @assert -Inf < a < b < Inf && abs(a) ≤ abs(b)
    @assert a ≤ 0 ≤ b || 0 < a < b

    Δ = (b - a) * middle(a, b)
    m = (1/√2) * expm1(-Δ) / (dawson(a/√2) * exp(-Δ) - dawson(b/√2))
    return clamp(m, a, b)
end

"""
    tnmom1i(a, b, μ, σ)

Mean of the normal distribution with variance -σ^2 and mean μ,
truncated to [a,b].
"""
tnmom1i(a, b, μ, σ) = μ + σ * tnmom1i((a - μ)/σ, (b - μ)/σ)
