using SpecialFunctions


function _F1(x::Real, y::Real; thresh=1e-7)
    @assert 0 < thresh < Inf
    #-Inf < x < Inf && -Inf < y < Inf || throw(DomainError())
    ϵ = exp(x^2 - y^2)
    if abs(x) > abs(y)
        _F1(y,x)
    elseif abs(x - y) ≤ thresh
        δ = y - x
        √π*x + (√π/2 + (-√π*x/6 + (-√π/12 + x*(√π/90 + (√π*x^2)/90)δ)δ)δ)δ
    elseif max(x,y) < 0
        (1 - ϵ) / (ϵ * erfcx(-y) - erfcx(-x))
    elseif min(x,y) > 0
        (1 - ϵ) / (erfcx(x) - ϵ * erfcx(y))
    else
        exp(-x^2) * (1 - ϵ) / (erf(y) - erf(x))
    end
end


function _F2(x::Real, y::Real; thresh=1e-7)
    @assert 0 < thresh < Inf
    -Inf < x < Inf && -Inf < y < Inf || throw(DomainError())
    ϵ = exp(x^2 - y^2)
    if abs(x) > abs(y)
        _F2(y,x)
    elseif abs(x - y) ≤ thresh
        δ = y - x
        √π*x^2 - √π/2 + (√π*x + (√π/3 - √π*x^2/3 + (((√π/30 + √π*x^2/45)x^2 - 4*√π/45)δ - √π*x/3)δ)δ)δ
    elseif max(x,y) < 0
        (x - ϵ * y) / (ϵ * erfcx(-y) - erfcx(-x))
    elseif min(x,y) > 0
        (x - ϵ * y) / (erfcx(x) - ϵ * erfcx(y))
    else
        exp(-x^2) * (x - ϵ * y) / (erf(y) - erf(x))
    end
end