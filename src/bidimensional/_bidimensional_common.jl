function _ellipse_quadratic(ρ::Real, x::Tuple{Real,Real})::Float64
    return x[1]^2 + x[2]^2 - 2ρ * x[1] * x[2]
end

function _log_normal_pdf(ρ::Real, x::Tuple{Real,Real})::Float64
    @assert -1 < ρ < 1
    #return -_ellipse_quadratic(ρ, x) / (2 * (1 - ρ^2)) - log(2π * sqrt(1 - ρ^2))
    return -_ellipse_quadratic(ρ, x) / (2 * (1 - ρ) * (1 + ρ)) - log(2π) - 0.5log(1 - ρ) - 0.5log(1 + ρ)
end