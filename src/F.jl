using Statistics, SpecialFunctions

"""
    _F1(x, y) = (exp(-y^2) - exp(-x^2)) / (erf(y) - erf(x))

without catastrophic cancellation. Defined only for x < y. 
_F1(-∞,+∞) is defined by taking the limit of the second argument first.
"""
function _F1(x::Real, y::Real)
    if !(x < y)
        throw(ArgumentError("Defined only for x < y; got x = $x, y = $y"))
    elseif abs(x) == abs(y)
        return zero(middle(x, y))
    elseif abs(x) > abs(y)
        return -_F1(-y, -x)
    end
    @assert x < y && abs(x) < abs(y)
    
    if isinf(y) # y == Inf
        return -inv(erfcx(x))
    elseif isinf(x) # x == -Inf
        inv(erfcx(-y))
    end
    
    Δm1 = expm1(x^2 - y^2)
    Δ = one(Δm1) + Δm1

    if 0 < x < y
        Δm1 / (erfcx(x) - Δ * erfcx(y))
    elseif x ≤ 0 ≤ y
        Δm1 * exp(-x^2) / (erf(y) - erf(x))
    end
end

"""
    _F4(x, y) = ((y+2x)exp(-y^2) - 3x exp(-x^2)) / (erf(y) - erf(x))

without catastrophic cancellation. _F1(±∞,±∞) is defined
by taking the limit of the second argument first.
"""
function _F4(x::Real, y::Real)  
    @assert 0 ≤ x ≤ y
    if isinf(y) && y > 0
        return -3x / erfcx(x)
    end
    
    Δm1 = expm1(x^2 - y^2)
    Δ = one(Δm1) + Δm1

    #return Δm1 * exp(-x^2) / (erf(y) - erf(x))

    #return Δm1 / (erfcx(x) - Δ * erfcx(y))

    return ((y+2x) * Δ - 3x) / (erfcx(x) - Δ * erfcx(y))
end

"""
    _F2(x, y) = (y * exp(-y^2) - x * exp(-x^2)) / (erf(y) - erf(x))

without catastrophic cancellation. _F2(±∞,±∞) is defined
by taking the limit in the second argument first.
"""
function _F2(x::Real, y::Real)
    if x > y && abs(x) > abs(y)
        return _F2(y, x)
    elseif x > y && abs(x) ≤ abs(y)
        return _F2(-x, -y)
    elseif x < y && abs(x) > abs(y)
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
without catastrophic cancellation. Defined only for x < y.
"""
function _F3(x::Real, y::Real)
    if !(x < y)
        throw(ArgumentError("Defined only for x < y; got x = $x, y = $y"))
    elseif abs(x) > abs(y)
        return _F3(-y, -x)
    end

    @assert x < y && abs(x) ≤ abs(y)

    if isinf(y)
        return -inv(erfcx(x))
    end

    Δ = exp(x^2 - y^2)

    if x < y < 0
        (Δ + 1) / (erfcx(-x) - Δ * erfcx(-y))
    elseif 0 < x < y
        (Δ + 1) / (Δ * erfcx(y) - erfcx(x))
    else
        (Δ + 1) * exp(-x^2) / (erf(x) - erf(y))
    end
end

"""
    _G2(x, y)

Computes (x * exp(-x^2) - y * exp(-y^2)) / (exp(-x^2) - exp(-y^2))
without catastrophic cancellation.
"""
function _G2(x::Real, y::Real)
    if x > y && abs(x) > abs(y)
        return _G2(y, x)
    elseif x > y && abs(x) ≤ abs(y)
        return -_G2(-x, -y)
    elseif x < y && abs(x) > abs(y)
        return -_G2(-y, -x)
    end

    @assert x ≤ y && abs(x) ≤ abs(y)

    if x == y
        return x - 1 / (2x)
    else
        Δm1 = expm1(x^2 - y^2)
        Δ = one(Δm1) + Δm1
        return (y * Δ - x) / Δm1
    end
end

"""
    _G3(x, y)

Computes (exp(-x^2) + exp(-y^2)) / (exp(-x^2) - exp(-y^2)).
"""
function _G3(x::Real, y::Real)
    if x > y && abs(x) > abs(y)
        return -_G3(y, x)
    elseif x > y && abs(x) ≤ abs(y)
        return _G3(-x, -y)
    elseif x < y && abs(x) > abs(y)
        return -_G3(-y, -x)
    end

    @assert x < y && abs(x) ≤ abs(y)

    Δm1 = expm1(x^2 - y^2)
    Δ = one(Δm1) + Δm1

    #(1 + Δ) / (1 - exp(x^2-y^2))

    return -(1 + Δ) / Δm1
end