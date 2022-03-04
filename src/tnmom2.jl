using Statistics, SpecialFunctions

"""
    tnmom2(a, b)

Second moment of the truncated standard normal distribution.
"""
function tnmom2(a::Real, b::Real)
    #return tnmom2c(0, a, b)

    if !(a ≤ b)
        return oftype(middle(a, b), NaN)
    elseif a == b
        return middle(a, b)^2
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
        ea = √(π/2) * erf(a / √2)
        eb = √(π/2) * erf(b / √2)
        fa = ea - a * exp(-a^2 / 2)
        fb = eb - b * exp(-b^2 / 2)
        m2 = (fb - fa) / (eb - ea)
        fb ≥ fa && eb ≥ ea || error("error: a=$a, b=$b")
        0 ≤ m2 ≤ 1 || error("error: a=$a, b=$b")
        return m2
    else # 0 ≤ a ≤ b
        exΔ = exp((a - b)middle(a, b))
        ea = √(π/2) * erfcx(a / √2)
        eb = √(π/2) * erfcx(b / √2)
        fa = ea + a
        fb = eb + b
        m2 = (fa - fb * exΔ) / (ea - eb * exΔ)
        a^2 ≤ m2 ≤ b^2 || error("error: a=$a, b=$b")
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
        #return σ^2 * tnmom2c(-μ / σ, α, β)
        return μ^2 + σ^2 * tnmom2(α, β) + 2μ * σ * tnmean(α, β)
    elseif iszero(σ) && a ≤ b
        # point mass
        # ⟹ if μ ∈ [a,b], 2nd moment is μ^2
        # ⟹ if μ < a, 2nd moment is a^2
        # ⟹ if μ > b, 2nd moment is b^2
        return clamp(μ / one(μ), a, b)^2
    else
        return oftype(middle(a, b), NaN)
    end
end

"""
    tnmom2i(a, b)

Second moment of the normal distribution with variance -1 and mean 0,
truncated to [a,b].
"""
function tnmom2i(a::Real, b::Real)
    if !(-Inf < a ≤ b < Inf)
        return oftype(middle(a, b), NaN)
    elseif a == b
        return middle(a, b)^2
    elseif abs(a) > abs(b)
        return tnmom2i(-b, -a)
    end

    @assert -Inf < a < b < Inf && abs(a) ≤ abs(b)
    @assert a ≤ 0 ≤ b || 0 < a < b

    Δ = (b - a) * middle(a, b)
    exΔ = exp(-Δ)
    da = dawson(a/√2)
    db = dawson(b/√2)
    m = ((da - a/√2)exΔ - (db - b/√2)) / (db - da * exΔ)
    return m
end

"""
    tnmom2i(a, b, μ, σ)

Second moment of the normal distribution with variance -σ^2 and mean μ,
truncated to [a,b].
"""
function tnmom2i(a, b, μ, σ)
    α = (a - μ) / σ
    β = (b - μ) / σ
    return μ^2 + σ^2 * tnmom2i(α, β) + 2μ * σ * tnmom1i(α, β)
end
