function _ellipse_quadratic(ρ::Real, x::Tuple{Real,Real})::Float64
    x[1]^2 + x[2]^2 - 2ρ * x[1] * x[2]
end

function _log_normal_pdf(ρ::Real, x::Tuple{Real,Real})::Float64
    @assert -1 < ρ < 1
    -_ellipse_quadratic(ρ, x) / (2 * (1 - ρ^2)) - log(2π * sqrt(1 - ρ^2))
end