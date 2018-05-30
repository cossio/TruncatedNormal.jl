import Cubature


"""
Returns val,C, such that the value of the integral is val * exp(C).
The constant C does not depend on f. 
"""
function _integrate_gauss_2D(f, ρ::Real, a::Tuple{Real,Real}, b::Tuple{Real,Real})::Tuple{Float64,Float64}
    @assert -1 < ρ < 1
    @assert -Inf < a[1] ≤ b[1] < Inf && -Inf < a[2] ≤ b[2] < Inf

    xm = _ellipse_rect_min_2D(ρ, a, b)
    C = _log_normal_pdf(ρ, xm)

    integrand(x) = f(x) * exp(_log_normal_pdf(ρ, (x[1], x[2])) - C)
    (val, err) = Cubature.hcubature(integrand, collect(a), collect(b))
    return val, C
end
