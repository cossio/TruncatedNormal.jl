import HCubature


"""
Integral of exp(H - |x|^2) * f(x) within the rectangle a .≤ x .≤ b.
"""
function _integrate(f, a::Vector, b::Vector; H::Real = 0, recursion::Integer = 3, Δ::Real=2.)
   _integrate(f, a, b, zeros(a), eye(eltype(a), length(a)); H = H, recursion = recursion, Δ = Δ)
end


"""
Integral of exp(H - (x-μ)ᵀΣinv(x-μ)) * f(x) within the rectangle a .≤ x .≤ b.
"""
function _integrate(f, a::Vector, b::Vector, μ::Vector, Σinv::Matrix; H::Real = 0, recursion::Integer = 3, Δ::Real=2.)
    @assert length(a) == length(b)
    @assert all(-Inf .< a .≤ b .< Inf)
    @assert recursion ≥ 0
    @assert 0 ≤ Δ < Inf
    @assert -Inf < H < Inf
    
    any(a .== b) && return 0., 0.
    
    # point closest to μ within the rectangle a .≤ x .≤ b.
    x0 = clamp.(μ, a, b)

    a0 = max.(a, x0 - Δ)
    b0 = min.(b, x0 + Δ)
    result, err = HCubature.hcubature(x -> exp(H - (x-μ)' * Σinv * (x-μ)) * f(x), a0, b0)

    a1 = copy(a0)
    b1 = copy(b0)
    if recursion > 0
        recursion -= 1
        for region = 1 : 3^length(a) - 1
            c = base(3, region, length(a))
            for i = 1 : length(a)
                if c[i] == '1'
                    a1[i] = a[i]
                    b1[i] = a0[i]
                elseif c[i] == '2'
                    a1[i] = b0[i]
                    b1[i] = b[i]
                else
                    a1[i] = a0[i]
                    b1[i] = b0[i]
                end
            end
            result1, err1 = _integrate(f, a1, b1, μ, Σinv; H=H, recursion=recursion, Δ=Δ)
            result += result1
            err += err1
        end
    end

    return result, err
end


"""
Integral of exp(H - (x-μ)ᵀΣinv(x-μ)) within the rectangle a .≤ x .≤ b.
"""
function _integrate(a::Vector, b::Vector, μ::Vector, Σinv::Matrix; H::Real = 0, recursion::Integer = 3, Δ::Real=2.)
    _integrate(x -> 1, a, b, μ, Σinv; H = H, recursion = recursion, Δ = 2)
end


"""
Integral of exp(H - |x|^2) within the rectangle a .≤ x .≤ b.
"""
function _integrate(a::Vector, b::Vector; H::Real = 0, recursion::Integer = 3, Δ::Real=2.)
    _integrate(x -> 1, a, b; H = H, recursion = recursion, Δ = Δ)
end


