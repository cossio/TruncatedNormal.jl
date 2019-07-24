using Statistics, SpecialFunctions
export tnmom1c

"""
    tnmom1c(c, a, b)

Computes ⟨(x - c)⟩ for the standard normal distribution truncated to the
interval [a, b]. This is more numerically stable than computing the mean
of x and subtracting c. Returns NaN if a > b.
"""
function tnmom1c(c::Real, a::Real, b::Real)
    if !(a ≤ b) || isnan(c)
        return oftype(middle(a, b), NaN)
    elseif a == b || isinf(c)
        return middle(a, b) - c
    elseif abs(a) > abs(b)
        return -tnmom1c(-c, -b, -a)
    elseif c ≤ 0
        #= At this point tnmom1(a,b) ≥ 0.
        Therefore the subtraction is fine. =#
        return tnmom1(a, b) - c
    elseif isinf(a) && isinf(b)
        return oftype(middle(a, b), -c)
    # elseif isfinite(a) && isinf(b)
    #     return √(2/π) / erfcx(a /√2) - c
    end

    @assert a < b
    @assert a ≤ 0 ≤ b || 0 ≤ a ≤ b
    @assert c > 0

    if a ≤ 0 ≤ b
        era = erf(a / √2)
        erb = erf(b / √2)
        exa = √(2/π) * exp(-a^2 / 2)
        exb = √(2/π) * exp(-b^2 / 2)
        if exa ≈ exb
            return -c
            fa = c * era
            fb = c * erb
        elseif era ≈ erb
            fa = exa
            fb = exb
        else
            fa = exa + c * era
            fb = exb + c * erb
        end
        return (fa - fb) / (erb - era)
    elseif 0 < a < b
        era = erfcx(a / √2)
        erb = erfcx(b / √2)
        if xerfcx_asym_thresh(a / √2, 3)
            fa = (1 - c/a) * √(2/π) + √2 * c/a * xerfcx_pi(a / √2)
        else
            fa = √(2/π) - c * era
        end
        if xerfcx_asym_thresh(b / √2, 3)
            fb = (1 - c/b) * √(2/π) + √2 * c/b * xerfcx_pi(b / √2)
        else
            fb = √(2/π) - c * erb
        end
        exΔ = exp(middle(a, b) * (a - b))
        return (fa - fb * exΔ) / (era - erb * exΔ)
    end
end

"""
    tnmom1c(c, a, b, μ, σ)

Computes ⟨(x - c)⟩ for the normal distribution with mean μ and standard
deviation σ, truncated to the interval [a, b]. This is more numerically stable
than computing the mean of x and subtracting c.
"""
function tnmom1c(c, a, b, μ, σ)
    α = (a - μ) / σ
    β = (b - μ) / σ
    γ = (c - μ) / σ
    return σ * tnmom1c(γ, α, β)
end
