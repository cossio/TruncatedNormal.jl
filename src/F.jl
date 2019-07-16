using Statistics, SpecialFunctions


"""
    _F1(x, y) = (exp(-y^2) - exp(-x^2)) / (erf(y) - erf(x))

without catastrophic cancellation. _F1(±∞,±∞) is defined
by taking the limit of the second argument first.
"""
function _F1(x::Real, y::Real)
    if x > y && abs(x) > abs(y)
        return _F1(y, x)
    elseif x > y && abs(x) ≤ abs(y)
        return -_F1(-x, -y)
    elseif x ≤ y && abs(x) > abs(y)
        return -_F1(-y, -x)
    end
    
    @assert x ≤ y && abs(x) ≤ abs(y)

    if x == y
        return -√π * x
    elseif isinf(y)
        return -sign(y) / erfcx(sign(y) * x)
    end
    
    Δm1 = expm1((x - y) * (x + y))
    Δ = one(Δm1) + Δm1

    if max(x, y) < 0
        Δm1 / (Δ * erfcx(-y) - erfcx(-x))
    elseif min(x, y) > 0 || y == Inf
        Δm1 / (erfcx(x) - Δ * erfcx(y))
    else
        Δm1 * exp(-x^2) / (erf(y) - erf(x))
    end
end

"""
    _F2(x, y) = (y * exp(-y^2) - x * exp(-x^2)) / (erf(y) - erf(x))

without catastrophic cancellation. _F2(±∞,±∞) is defined
by taking the limit in the second argument first.
"""
function _F2(x::Real, y::Real)::Float64
    if x > y && abs(x) > abs(y)
        return _F2(y, x)
    elseif x > y && abs(x) ≤ abs(y)
        return _F2(-x, -y)
    elseif x ≤ y && abs(x) > abs(y)
        return _F2(-y, -x)
    end

    @assert x ≤ y && abs(x) ≤ abs(y)

    if x == -Inf && y == Inf
        return zero(x + y)
    elseif isinf(y)
        return -sign(y) * x / erfcx(sign(y) * x)
    elseif x == y
        return √π * (1 - 2x^2) / 2
    end

    Δm1 = expm1(x^2 - y^2)
    Δ = Δm1 + 1

    if max(x, y) < 0
        (Δ * y - x) / (Δ * erfcx(-y) - erfcx(-x))
    elseif min(x, y) > 0
        (Δ * y - x) / (erfcx(x) - Δ * erfcx(y))
    else
        exp(-x^2) * (Δ * y - x) / (erf(y) - erf(x))
    end
end

"""
    _F3(x, y) = 
    
Computes (exp(-x^2) + exp(-y^2)) / (erf(x) - erf(y))
without catastrophic cancellation.
"""
function _F3(x::Real, y::Real)::Float64
    if x > y && abs(x) > abs(y)
        return -_F3(y, x)
    elseif x > y && abs(x) ≤ abs(y)
        return -_F3(-x, -y)
    elseif x < y && abs(x) > abs(y)
        return _F3(-y, -x)
    end

    @assert x < y && abs(x) ≤ abs(y)

    if x == -Inf && y == Inf
        return zero(x + y)
    elseif isinf(y)
        return -inv(erfcx(x))
    end

    Δ = exp(x^2 - y^2)

    if max(x, y) < 0
        (Δ + 1) / (erfcx(-x) - Δ * erfcx(-y))
    elseif min(x, y) > 0
        (Δ + 1) / (Δ * erfcx(y) - erfcx(x))
    else
        (Δ + 1) * exp(-x^2) / (erf(x) - erf(y))
    end
end
