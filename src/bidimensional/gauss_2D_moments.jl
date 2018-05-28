"⟨x^i * x^j⟩"
function _gauss2Dmoment(ρ::Real, a::Tuple{Real,Real}, b::Tuple{Real,Real}, i::Real, j::Real)
    _gauss2Dmoment_nonorm(ρ, a, b, i, j) / _gauss2Dtrunc_partition(ρ, a, b)
end

function _gauss2Dtrunc_partition(ρ::Real, a::Tuple{Real,Real}, b::Tuple{Real,Real})
    _gauss2Dmoment_nonorm(ρ, a, b, 0, 0)
end

function _gauss2Dmoment_nonorm(ρ::Real, a::Tuple{Real,Real}, b::Tuple{Real,Real}, i::Real, j::Real)
    _integrate_gauss_2D(x -> x[1]^i * x[2]^j, ρ, a, b)
end

"⟨x^i * x^j⟩"
function gauss2Dmoment(μ::Tuple{Real,Real}, Σ::AbstractMatrix{T}, a::Tuple{Real,Real}, b::Tuple{Real,Real}, m::Real, n::Real) where {T <: Real}

    @assert size(Σ) == (2,2) && Σ[1,2] == Σ[2,1]

    α = ((a[1] - μ[1]) / √Σ[1,1], (a[2] - μ[2]) / √Σ[2,2])
    β = ((b[1] - μ[1]) / √Σ[1,1], (b[2] - μ[2]) / √Σ[2,2])
    ρ = Σ[1,2] / √(Σ[1,1] * Σ[2,2])
    @assert -1 < ρ < 1

    Z = _gauss2Dtrunc_partition(ρ, α, β)

    sum(binomial(m,i) * binomial(n,j) * μ[1]^(m-i) * μ[2]^(n-j) * √(Σ[1,1]^i * Σ[2,2]^j) * 
        _gauss2Dmoment_nonorm(ρ, α, β, i, j) / Z
        for i = 0 : m for j = 0 : n)
end

"⟨x1⟩"
function gauss2Dmoment1(μ::Tuple{Real,Real}, Σ::AbstractMatrix{T}, a::Tuple{Real,Real}, b::Tuple{Real,Real}) where {T <: Real}
    gauss2Dmoment(μ, Σ, a, b, 1, 0)
end

"⟨x2⟩"
function gauss2Dmoment2(μ::Tuple{Real,Real}, Σ::AbstractMatrix{T}, a::Tuple{Real,Real}, b::Tuple{Real,Real}) where {T <: Real}
    gauss2Dmoment(μ, Σ, a, b, 0, 1)
end

"⟨x1^2⟩"
function gauss2Dmoment11(μ::Tuple{Real,Real}, Σ::AbstractMatrix{T}, a::Tuple{Real,Real}, b::Tuple{Real,Real}) where {T <: Real}
    gauss2Dmoment(μ, Σ, a, b, 2, 0)
end

"⟨x2^2⟩"
function gauss2Dmoment22(μ::Tuple{Real,Real}, Σ::AbstractMatrix{T}, a::Tuple{Real,Real}, b::Tuple{Real,Real}) where {T <: Real}
    gauss2Dmoment(μ, Σ, a, b, 0, 2)
end

"⟨x1*x2⟩"
function gauss2Dmoment12(μ::Tuple{Real,Real}, Σ::AbstractMatrix{T}, a::Tuple{Real,Real}, b::Tuple{Real,Real}) where {T <: Real}
    gauss2Dmoment(μ, Σ, a, b, 1, 1)
end

"Returns the truncated μ (means) and Σ (covariance matrix)"
function gauss2Dtruncstats(μ::Tuple{Real,Real}, Σ::AbstractMatrix{T}, a::Tuple{Real,Real}, b::Tuple{Real,Real}) where {T <: Real}

    @assert size(Σ) == (2,2) && Σ[1,2] == Σ[2,1]

    α = ((a[1] - μ[1]) / √Σ[1,1], (a[2] - μ[2]) / √Σ[2,2])
    β = ((b[1] - μ[1]) / √Σ[1,1], (b[2] - μ[2]) / √Σ[2,2])
    ρ = Σ[1,2] / √(Σ[1,1] * Σ[2,2])
    @assert -1 < ρ < 1

    Z = _gauss2Dtrunc_partition(ρ, α, β)

    ξ1 = _gauss2Dmoment_nonorm(ρ, α, β, 1, 0) / Z
    ξ2 = _gauss2Dmoment_nonorm(ρ, α, β, 0, 1) / Z
    ξ11 = _gauss2Dmoment_nonorm(ρ, α, β, 2, 0) / Z
    ξ22 = _gauss2Dmoment_nonorm(ρ, α, β, 0, 2) / Z
    ξ12 = _gauss2Dmoment_nonorm(ρ, α, β, 1, 1) / Z

    tμ1 = μ[1] + √Σ[1,1] * ξ1
    tμ2 = μ[2] + √Σ[2,2] * ξ2
    tΣ11 = μ[1]^2 + 2μ[1] * √Σ[1,1] * ξ1 + Σ[1,1] * ξ11
    tΣ22 = μ[2]^2 + 2μ[2] * √Σ[2,2] * ξ2 + Σ[2,2] * ξ22
    tΣ12 = μ[1]μ[2] + μ[1] * √Σ[2,2] * ξ2 + μ[2] * √Σ[1,1] * ξ1 + √(Σ[1,1]Σ[2,2]) * ξ12
    
    return (tμ1,tμ2), [tΣ11 tΣ12; tΣ12 tΣ22]
end


function gauss2Dtruncstats(μ::Vector{T}, Σ::AbstractMatrix{T}, a::Vector{T}, b::Vector{T}) where {T <: Real}
    @assert length(μ) == length(a) == length(b) == 2
    @assert size(Σ) == (2,2)
    tμ, tΣ = gauss2Dtruncstats((μ[1],μ[2]), Σ, (a[1],a[2]), (b[1],b[2]))
    return [tμ[1], tμ[2]], tΣ
end