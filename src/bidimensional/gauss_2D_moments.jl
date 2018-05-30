"⟨x^i * x^j⟩"
function _gauss2Dmoment(ρ::Real, a::Tuple{Real,Real}, b::Tuple{Real,Real}, i::Real, j::Real)
    v1, C1 = _gauss2Dmoment_nonorm(ρ, a, b, i, j)
    v0, C0 = _gauss2Dtrunc_partition(ρ, a, b)
    @assert C0 == C1
    v1 / v0
end

function _gauss2Dtrunc_partition(ρ::Real, a::Tuple{Real,Real}, b::Tuple{Real,Real})
    _gauss2Dmoment_nonorm(ρ, a, b, 0, 0)
end

function _gauss2Dmoment_nonorm(ρ::Real, a::Tuple{Real,Real}, b::Tuple{Real,Real}, i::Real, j::Real)
    _integrate_gauss_2D(x -> x[1]^i * x[2]^j, ρ, a, b)
end

"⟨x1^m * x2^n⟩"
function gauss2Dmoment(μ::Tuple{Real,Real}, Σ::AbstractMatrix{T}, a::Tuple{Real,Real}, b::Tuple{Real,Real}, m::Integer, n::Integer) where {T <: Real}

    @assert size(Σ) == (2,2) && Σ[1,2] == Σ[2,1]

    α = ((a[1] - μ[1]) / √Σ[1,1], (a[2] - μ[2]) / √Σ[2,2])
    β = ((b[1] - μ[1]) / √Σ[1,1], (b[2] - μ[2]) / √Σ[2,2])
    ρ = Σ[1,2] / √(Σ[1,1] * Σ[2,2])
    @assert -1 < ρ < 1

    # The limits ρ → ± 1 can be computed analytically..... but this turned out to be worse
    # than just letting the integral do its job
    # if ρ ≤ -0.99
    #     if α[1] ≤ α[2]
    #         ξ = (i,j) -> α[1]^i * α[2]^j
    #     elseif β[1] ≥ β[2]
    #         ξ = (i,j) -> β[1]^i * β[2]^j
    #     else
    #         ξ = (i,j) -> (-1)^j / (i + j + 1) * (min(β[1], -α[2])^(i + j + 1) - max(α[1], -β[2])^(i + j + 1)) / (min(β[1], -α[2]) - max(α[1], -β[2]))
    #     end
    # elseif ρ ≥ 0.99
    #     if α[1] ≤ β[2]
    #         ξ = (i,j) -> α[1]^i * β[2]^j
    #     elseif β[1] ≤ α[2]
    #         ξ = (i,j) -> β[1]^i * α[2]^j
    #     else
    #         ξ = (i,j) -> 1 / (i + j + 1) * (min(β[1], β[2])^(i + j + 1) - max(α[1], α[2])^(i + j + 1)) / (min(α[1], α[2]) - max(β[1], β[2]))
    #     end
    # else
    #     Z = _gauss2Dtrunc_partition(ρ, α, β)
    #     ξ = (i,j) -> _gauss2Dmoment_nonorm(ρ, α, β, i, j) / Z
    # end

    Z, C0 = _gauss2Dtrunc_partition(ρ, α, β)

    ξ = (i,j) -> begin
        v1, C1 = _gauss2Dmoment_nonorm(ρ, α, β, i, j)
        v1 / Z
    end 

    sum(binomial(m,i) * binomial(n,j) * μ[1]^(m-i) * μ[2]^(n-j) * √(Σ[1,1]^i * Σ[2,2]^j) * ξ(i,j)
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

    Z, _ = _gauss2Dtrunc_partition(ρ, α, β)

    ξ10, _ = _gauss2Dmoment_nonorm(ρ, α, β, 1, 0)
    ξ01, _ = _gauss2Dmoment_nonorm(ρ, α, β, 0, 1)
    ξ20, _ = _gauss2Dmoment_nonorm(ρ, α, β, 2, 0)
    ξ02, _ = _gauss2Dmoment_nonorm(ρ, α, β, 0, 2)
    ξ11, _ = _gauss2Dmoment_nonorm(ρ, α, β, 1, 1)

    ξ10 /= Z
    ξ01 /= Z
    ξ20 /= Z
    ξ02 /= Z
    ξ11 /= Z

    tμ1 = μ[1] + √Σ[1,1] * ξ10
    tμ2 = μ[2] + √Σ[2,2] * ξ01
    tΣ11 = μ[1]^2 + 2μ[1] * √Σ[1,1] * ξ10 + Σ[1,1] * ξ20
    tΣ22 = μ[2]^2 + 2μ[2] * √Σ[2,2] * ξ01 + Σ[2,2] * ξ02
    tΣ12 = μ[1]μ[2] + μ[1] * √Σ[2,2] * ξ01 + μ[2] * √Σ[1,1] * ξ10 + √(Σ[1,1]Σ[2,2]) * ξ11
    
    return (tμ1,tμ2), [tΣ11 tΣ12; tΣ12 tΣ22]
end


function gauss2Dtruncstats(μ::Vector{T}, Σ::AbstractMatrix{T}, a::Vector{T}, b::Vector{T}) where {T <: Real}
    @assert length(μ) == length(a) == length(b) == 2
    @assert size(Σ) == (2,2)
    tμ, tΣ = gauss2Dtruncstats((μ[1],μ[2]), Σ, (a[1],a[2]), (b[1],b[2]))
    return [tμ[1], tμ[2]], tΣ
end