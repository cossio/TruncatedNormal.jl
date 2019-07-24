using Statistics, SpecialFunctions
export tnmom2c

"""
    tnmom2c(c, a, b)

Computes ⟨(x - c)^2⟩ for the standard normal distribution truncated to the
interval [a, b].
"""
function tnmom2c(c::Real, a::Real, b::Real)
    if !(a ≤ b) || isnan(c)
        return oftype(middle(a, b), NaN)
    elseif a == b
        return oftype(middle(a, b), (a - c) * (b - c))
    elseif isinf(c)
        return c^2
    elseif abs(a) > abs(b)
        return tnmom2c(-c, -b, -a)
    elseif isinf(a) && isinf(b)
        return oftype(middle(a, b), 1 + c^2)
    elseif c ≤ 0
        #= at this point tnmom1(a, b) ≥ 0, so the subtraction is fine =#
        m1 = tnmom1(a, b)
        @assert m1 ≥ 0
        return tnmom2(a, b) + c^2 - 2c * m1
    elseif isfinite(a) && isinf(b)
        return 1 + c^2 + √(2/π) * (a - 2c) / erfcx(a / √2)
    end

    @assert a < b
    @assert 0 ≤ a ≤ b || a ≤ 0 ≤ b

    exp_a = exp(-a^2 / 2)
    exp_b = exp(-b^2 / 2)

    erfcx_abs_a = erfcx(abs(a) / √2)
    erfcx_abs_b = erfcx(abs(b) / √2)
    erfc_abs_a = exp_a * erfcx_abs_a
    erfc_abs_b = exp_b * erfcx_abs_b

    ra = √(2/π) * (a - 2c) / (1 + c^2)    # a / (1 + c^2) - 2 / (1/c + c)
    rb = √(2/π) * (b - 2c) / (1 + c^2)

    Δ = (b - a) * middle(a, b)

    if a ≤ 0 ≤ b
        erfc_a = 2 - erfc_abs_a
        erfc_b = erfc_abs_b
        erfc_ma = erfc_abs_a
        erfc_mb = 2 - erfc_abs_b
        erf_a = erfc_abs_a - 1
        erf_b = 1 - erfc_abs_b

        γa = ra * exp_a - erf_a
        γb = rb * exp_b - erf_b
        @assert γa ≥ γb

        if γa ≥ 0 ≥ γb
            m = γa - γb
        elseif γa ≤ 0 && γb ≤ 0
            m = (ra * exp_a + erfc_a) - (rb * exp_b + erfc_b)
        elseif γa ≥ 0 && γb ≥ 0
            m = (ra * exp_a - erfc_ma) - (rb * exp_b - erfc_mb)
        end
        m /= (erf_b - erf_a)
    elseif 0 ≤ a ≤ b
        erfc_a = erfc_abs_a
        erfc_b = erfc_abs_b
        erfc_ma = 2 - erfc_abs_a
        erfc_mb = 2 - erfc_abs_b
        erf_a = 1 - erfc_a
        erf_b = 1 - erfc_b
        erfcx_a = erfcx_abs_a
        erfcx_b = erfcx_abs_b

        γa = ra * exp_a - erf_a
        γb = rb * exp_b - erf_b
        @assert γa ≥ γb

        if γa ≥ 0 ≥ γb
            m = (γa - γb) / (erfc_a - erfc_b)
        elseif γa ≤ 0 && γb ≤ 0 # ζ ≤ a ≤ b
            exp_Δ = exp(-Δ)
            m = (ra + erfcx_a) - (rb + erfcx_b) * exp_Δ
            m /= erfcx_a - erfcx_b * exp_Δ
        elseif γa ≥ 0 && γb ≥ 0
            m = (ra * exp_a - erfc_ma) - (rb * exp_b - erfc_mb)
            m /= erfc_a - erfc_b
        else
            error("trouble")
        end
    else
        error("trouble")
    end

    return (1 + c^2) * m

    # γa = √(2/π) * (a - 2c) / (1 + c^2) * exp(-a^2 / 2) - ea
    # γb = √(2/π) * (b - 2c) / (1 + c^2) * exp(-b^2 / 2) - eb
    #
    # ea = erf(a / √2)
    # eb = erf(b / √2)
    #
    # if γb ≤ 0 ≤ γa
    #     return (1 + c^2) * (γa - γb) / (eb - ea)
    # else
    #     ea = erfcx(a / √2)
    #     eb = erfcx(b / √2)
    #     γa = √(2/pi) * (a - 2c) / (1 + c^2) + ea
    #     γb = √(2/pi) * (b - 2c) / (1 + c^2) + eb
    #     @assert 0 ≤ Δ ≤ 1
    #     #@assert γb ≤ γa
    #     @assert eb ≤ ea
    #     return (1 + c^2) * (γa - γb * Δ) / (ea - eb * Δ)
    # end
end

"""
    tnmom2c(c, a, b, μ, σ)

Computes ⟨(x - c)^2⟩ for the normal distribution with mean μ and standard
deviation σ, truncated to the interval [a, b].
"""
function tnmom2c(c, a, b, μ, σ)
    α = (a - μ) / σ
    β = (b - μ) / σ
    γ = (c - μ) / σ
    return σ^2 * tnmom2c(γ, α, β)
end
