using Statistics, SpecialFunctions
export tnmom1, tnmom2, tnvar, tnmean

"""
    tnmom2(a, b)

Second moment of the truncated standard normal distribution.
"""
function tnmom2(a::Real, b::Real)
    #return tnmom2c(0, a, b)

    if !(a ≤ b)
        return oftype(middle(a, b), NaN)
    elseif a == b
        return oftype(middle(a, b), a * b)
    elseif abs(a) > abs(b)
        return tnmom2(-b, -a)
    elseif isinf(a) && isinf(b)
        return one(middle(a, b))
    elseif isinf(b)
        return 1 + √(2 / π) * a / erfcx(a / √2)
    end

    @assert a < b < Inf && abs(a) ≤ abs(b)
    @assert a ≤ 0 ≤ b || 0 ≤ a ≤ b

    if a ≤ 0 ≤ b
        ea = erf(a / √2)
        eb = erf(b / √2)
        fa = ea - √(2/π) * a * exp(-a^2 / 2)
        fb = eb - √(2/π) * b * exp(-b^2 / 2)
        m2 = (fb - fa) / (eb - ea)
        @assert fb ≥ fa && eb ≥ ea
        @assert 0 ≤ m2 ≤ 1
        return m2
    else # 0 ≤ a ≤ b
        exΔ = exp((a - b)middle(a, b))
        ea = erfcx(a / √2)
        eb = erfcx(b / √2)
        fa = ea + √(2/π) * a
        fb = eb + √(2/π) * b
        m2 = (fa - fb * exΔ) / (ea - eb * exΔ)
        @assert a^2 ≤ m2 ≤ b^2
        return m2
    end
end

"""
    tnmom2(a, b, μ, σ)

Second moment of the truncated normal distribution, where μ, σ are the mean and
standard deviation of the untruncated distribution.
"""
function tnmom2(a, b, μ, σ)
    if σ > 0
        α = (a - μ) / σ
        β = (b - μ) / σ
        return tnmom2c(-μ / σ, α, β)
    elseif iszero(σ) && a ≤ b
        return clamp(μ^2 / one(μ), a, b)
    else
        return oftype(middle(a, b), NaN)
    end
end
