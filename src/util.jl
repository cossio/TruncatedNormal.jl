using Statistics, SpecialFunctions

"""
    xerfcx_pi(x)

Computes 1/√π - x * erfcx(x), without cancellation.
"""
function xerfcx_pi(x::Real)
    if xerfcx_asym_thresh(x, 3)
        return xerfcx_asym_pi(x, 3)
    else
        return 1/√π - x * erfcx(x)
    end
end

"""
    xerfcx_asym(x, m)

Evaluates x * erfcx(x) by its asymptotic series for large x > 0,
up to m'th order.
"""
function xerfcx_asym(x::Real, m::Integer)
    #= NIST Handbook of Mathematical Functions (2010) 1st Edition, 7.12.1 =#
    @assert m ≥ 0 && x > 0
    s = one(x / 1)
    for j = m : -1 : 0
        #s = 1 - (j + 1/2) / x^2 * s
        s = muladd((j + 1/2) / x / x, -s, 1)
    end
    return s / √π
end

"""
    xerfcx_asym(x, m)

Estimated reminder of xerfcx_asym(x, m).
"""
function xerfcx_asym_rem(x::Real, m::Integer)
    #= NIST Handbook of Mathematical Functions (2010) 1st Edition, 7.12.1. =#
    @assert m ≥ 0 && x > 0
    return 1 / √π * exp(lgamma(m + 1/2) - (2m + 1) * x) / erfcx(x)
end

"""
    xerfcx_asym_pi(x, m)

Computes 1/√π - x * erfcx(x) by its asymptotic series for large x > 0, up to
m'th order.
"""
function xerfcx_asym_pi(x::Real, m::Integer)
    #= NIST Handbook of Mathematical Functions (2010) 1st Edition, 7.12.1 =#
    @assert m ≥ 0 && x > 0
    s = one(x / 1)
    for j = m : -1 : 1
        #s = 1 - s * (j + 1/2) / x^2
        s = muladd((j + 1/2) / x / x, -s, 1)
    end
    s = s / 2 / x / x
    return s / √π
end

"""
    xerfcx_asym_thresh(x, m)

Returns true if the asymptotic series of x * erfcx(x) to order m can be used
for this value of x without appreciable error.
"""
xerfcx_asym_thresh(x::Real, m::Integer) = xerfcx_asym_thresh(x, Val(m))
function xerfcx_asym_thresh(x::T, ::Val{m}) where {T <: Real, m}
    @assert m isa Integer
    @assert m > 0
    return x > (√π * eps(float(T)))^(-1/2m) * exp(lgamma(m + 1/2) / 2m)
end
xerfcx_asym_thresh(x::Union{Float64,Integer}, ::Val{3}) = x > 452
xerfcx_asym_thresh(x::Union{Float32}, ::Val{3}) = x > 16
xerfcx_asym_thresh(x::Union{Float16}, ::Val{3}) = x > 4
