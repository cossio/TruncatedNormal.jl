export erf, erfc, inverfc

#= Numerical recipes implementation of Error function
Press, et al (2007) Numerical Recipes: The Art of Scientific Computing. 3rd Ed.
Section 6.2.2 =#

const erfcof = (-1.3026537197817094, 6.4196979235649026e-1,
    1.9476473204185836e-2,-9.561514786808631e-3,-9.46595344482036e-4,
    3.66839497852761e-4,4.2523324806907e-5,-2.0278578112534e-5,
    -1.624290004647e-6,1.303655835580e-6,1.5626441722e-8,-8.5238095915e-8,
    6.529054439e-9,5.059343495e-9,-9.91364156e-10,-2.27365122e-10,
    9.6467911e-11, 2.394038e-12,-6.886027e-12,8.94487e-13, 3.13092e-13,
    -1.12708e-13,3.81e-16,7.106e-15,-1.523e-15,-9.4e-17,1.21e-16,-2.8e-17)

function erf(x::Float64)
    if x ≥ 0
        return 1 - erfccheb(x)
    else
        return erfccheb(-x) - 1
    end
end

erf(x::Float64, y::Float64) = erfc(y, x)

function erfc(x::Float64)
    if x ≥ 0
        return erfccheb(x)
    else
        return 2 - erfccheb(-x)
    end
end

"""
    erfc(y) - erfc(x)
"""
function erfc(x::Float64, y::Float64)
    if x > y
        return -erfc(y, x)
    elseif abs(x) > abs(y)
        return -erfc(-x, -y)
    elseif x < 0 ≤ y
        return erf(x) - erf(y)
    elseif 0 ≤ x ≤ y
        return erfccheb(x, y)
    else
        return NaN
    end
end

function erfcx(x::Float64)
    if x ≥ 0
        return erfcxcheb(x)
    else
        return 2exp(x^2) - erfcxcheb(-x)
    end
end

inverf(p::Float64) = inverfc(1 - p)

# erfc
function erfccheb(z::Float64)
    z ≥ 0 || throw(ArgumentError("erfccheb requires nonnegative argument"))
    d = dd = 0.0
    t = 2 / (2 + z)
    ty = 4t - 2
    for j = length(erfcof) : -1 : 2
        tmp = d
        d = ty * d - dd + erfcof[j]
        dd = tmp
    end
    return t * exp(-z^2 + 0.5 * (erfcof[1] + ty * d) - dd)
end

# erfcx
function erfcxcheb(z::Float64)
    if !(z ≥ 0)
        return NaN
    end
    d = dd = 0.0
    t = 2 / (2 + z)
    ty = 4t - 2
    for j = length(erfcof) : -1 : 2
        tmp = d
        d = ty * d - dd + erfcof[j]
        dd = tmp
    end
    return t * exp(0.5 * (erfcof[1] + ty * d) - dd)
end

# erfc(y) - erfc(x)
function erfccheb(x::Float64, y::Float64)
    if !(0 ≤ x ≤ y)
        return NaN
    end
    dx = ddx = 0.0
    dy = ddy = 0.0
    tx = 2 / (2 + x)
    ty = 2 / (2 + y)
    qx = 4tx - 2
    qy = 4ty - 2
    for j = length(erfcof) : -1 : 2
        tmpx = dx
        tmpy = dy
        dx = qx * dx - ddx + erfcof[j]
        dy = qy * dy - ddy + erfcof[j]
        ddx = tmpx
        ddy = tmpy
    end
    fx = tx * exp(-x^2 + (erfcof[1] + qx * dx) / 2 - ddx)
    fy = ty * exp(-y^2 + (erfcof[1] + qy * dy) / 2 - ddy)
    return fy - fx
end

function inverfc(p::Float64)
    if p ≥ 2.0
        return -100.0
    elseif p <= 0.0
        return 100.0
    end

    pp = p < 1.0 ? p : 2.0 - p
    t = sqrt(-2 * log(pp / 2))
    x = -0.70711 * ((2.30753 + 0.27061t) / (1 + (0.99229 + 0.04481t)t) - t)
    for j = 1:2
        err = erfc(x) - pp
        x += err/(1.12837916709551257exp(-x^2) - x * err)
    end
    return p < 1.0 ? x : -x
end
