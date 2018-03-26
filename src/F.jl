using SpecialFunctions


"""
_F1(x, y) =

    (exp(-x^2) - exp(-y^2)) / (erf(y) - erf(x))

without catastrophic cancellation. _F1(±∞,±∞) is defined
by taking the limit of the second argument first.
"""
function _F1(x::Real, y::Real; thresh=1e-7)  
    if abs(x) > abs(y)
        return _F1(y, x)
    elseif isinf(y)
        return sign(y) / erfcx(sign(y) * x)
    elseif abs(x - y) ≤ thresh
        ϵ = y - x
        return √π*x + (√π/2 + (-√π*x/6 + (-√π/12 + x*(√π/90 + (√π*x^2)/90)ϵ)ϵ)ϵ)ϵ
    end
    
    Δ = exp(x^2 - y^2)

    if max(x, y) < 0
        (1 - Δ) / (Δ * erfcx(-y) - erfcx(-x))
    elseif min(x, y) > 0 || y == Inf
        (1 - Δ) / (erfcx(x) - Δ * erfcx(y))
    else
        exp(-x^2) * (1 - Δ) / (erf(y) - erf(x))
    end
end


"""
_F2(x, y) =

    (x * exp(-x^2) - y * exp(-y^2)) / (erf(y) - erf(x))

without catastrophic cancellation. _F2(±∞,±∞) is defined
by taking the limit in the second argument first.
"""
function _F2(x::Real, y::Real; thresh=1e-7)::Float64
    if abs(x) > abs(y)
        return _F2(y, x)
    elseif x == Inf && y == -Inf || x == -Inf && y == Inf
        return 0.
    elseif isinf(y)
        return sign(y) * x / erfcx(sign(y) * x)
    elseif abs(x - y) ≤ thresh
        ϵ = y - x
        return √π*x^2 - √π/2 + (√π*x + (√π/3 - √π*x^2/3 + (((√π/30 + √π*x^2/45)x^2 - 4*√π/45)ϵ - √π*x/3)ϵ)ϵ)ϵ
    end

    Δ = exp(x^2 - y^2)
    
    if max(x, y) < 0
        (x - Δ * y) / (Δ * erfcx(-y) - erfcx(-x))
    elseif min(x, y) > 0
        (x - Δ * y) / (erfcx(x) - Δ * erfcx(y))
    else
        exp(-x^2) * (x - Δ * y) / (erf(y) - erf(x))
    end
end
