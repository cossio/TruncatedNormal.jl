"""
<f> under multivariate normal
"""
function mtnexpectation(f, a::Vector, b::Vector, μ::Vector, Σinv::Matrix)
    @assert length(a) == length(b) == length(μ)
    @assert size(Σinv) == (length(a), length(a))
    @assert all(-Inf .< a .≤ b .< Inf)  # I do not support infinite bounds yet ... maybe in the future
    @assert all(-Inf .< μ .< Inf)

    x0 = clamp.(μ, a, b)
    fmax = exp(-(x0 - μ)' * Σinv * (x0 - μ))

    num, errn = _integrate(f, a, b, μ, Σinv; H=fmax)
    den, errd = _integrate(   a, b, μ, Σinv; H=fmax)

    return num / den
end

mtnexpectation(f, a::Vector, b::Vector) = mtnexpectation(f, a, b, zeros(a), eye(eltype(a), length(a)))


"""
mean of the truncated multivariate normal
"""
mtnmean(a::Vector, b::Vector, μ::Vector, Σinv::Matrix) = mtnexpectation(identity, a, b, μ, Σinv)
mtnmean(a::Vector, b::Vector) = mtnmean(a, b, zeros(a), eye(eltype(a), length(a)))


"""
matrix <xi*xj> under truncated multivariate normal
"""
function mtnpairprod(a::Vector, b::Vector, μ::Vector, Σinv::Matrix)
    f(x) = [x[i] * x[j] for i = 1 : length(a), j = 1 : length(b)]
    return mtnexpectation(f, a, b, μ, Σinv)
end

mtnpairprod(a::Vector, b::Vector) = mtnpairprod(a, b, zeros(a), eye(eltype(a), length(a)))

"""
covariance matrix of truncated multivariate normal
"""
function mtncov(a::Vector, b::Vector, μ::Vector, Σinv::Matrix)
    m = mtnmean(a, b, μ, Σinv)
    f(x) = [(x[i] - m[i]) * (x[j] - m[j]) for i = 1 : length(a), j = 1 : length(b)]
    return mtnexpectation(f, a, b, μ, Σinv)
end

mtncov(a::Vector, b::Vector) = mtncov(a, b, zeros(a), eye(eltype(a), length(a)))