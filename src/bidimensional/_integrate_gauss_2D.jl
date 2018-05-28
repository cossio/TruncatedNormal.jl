import Cubature

function _integrate_gauss_2D(f, ρ::Real, a::Tuple{Real,Real}, b::Tuple{Real,Real}, C::Real)::Float64
    @assert -1 < ρ < 1
    @assert -Inf < a[1] ≤ b[1] < Inf && -Inf < a[2] ≤ b[2] < Inf
    integrand(x) = f(x) * exp(_log_normal_pdf(ρ, (x[1], x[2])) - C)
    (val, err) = Cubature.hcubature(integrand, collect(a), collect(b))
    return val
end


function _integrate_gauss_2D(f, ρ::Real, a::Tuple{Real,Real}, b::Tuple{Real,Real})::Float64
    xm = _ellipse_rect_min_2D(ρ, a, b)
    C = _log_normal_pdf(ρ, xm)
    return _integrate_gauss_2D(f, ρ, a, b, C)
end
