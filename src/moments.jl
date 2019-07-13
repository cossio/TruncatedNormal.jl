export tnmean, tnvar, tnmom1, tnmom2

"""
    tnmean(a,b)

Mean of the truncated standard normal distribution.
"""
function tnmean(a, b)
    if a < b
        √(2/π) * _F1(a / √2, b / √2)
    elseif a == b
        a
    else
        throw(ArgumentError("a must be < b; got a = $a, b = $b"))
    end
end

"""
    tnmean(a, b, μ, σ)

Mean of the truncated normal distribution, where μ, σ
are the mean and standard deviation of the untruncated
distribution.
"""
function tnmean(a, b, μ, σ)
    α = (a - μ) / σ
    β = (b - μ) / σ
    return μ + tnmean(α, β) * σ
end

"""
    tnmom1(a, b)

Mean of the truncated standard normal distribution.

    tnmom1(a, b, μ, σ)

Mean of the truncated normal distribution, where μ, σ
are the mean and standard deviation of the untruncated
distribution.
"""
function tnmom1 end
tnmom1 = tnmean

"""
    tnmom2(a, b)

Second moment of the truncated standard normal distribution.
"""
function tnmom2(a, b)
    if a < b
        1 + (2/√π) * _F2(a / √2, b / √2)
    elseif a == b
        a^2
    else
        throw(ArgumentError("a must be < b; got a = $a, b = $b"))
    end
end

"""
    tnmom2(a, b, μ, σ)

Second moment of the truncated normal distribution, where μ, σ
are the mean and standard deviation of the untruncated distribution.
"""
function tnmom2(a, b, μ, σ)
    α = (a - μ) / σ
    β = (b - μ) / σ
    μ * (μ + 2σ * tnmom1(α, β)) + tnmom2(α, β) * σ^2
end

"""
    tnvar(a,b)

Variance of the truncated standard normal distribution
"""
function tnvar(a, b)
    return 1 + 2/√π * _F2(a/√2, b/√2) - 2/π * _F1(a/√2, b/√2)^2
    if a < b
        1 + 2/√π * _F2(a/√2, b/√2) - 2/π * _F1(a/√2, b/√2)^2
    elseif a == b
        a
    else
        throw(ArgumentError("a must be ≤ b"))
    end
end

"""
    tnvar(a, b, μ, σ)

Variance of the truncated normal distribution, where μ, σ are the mean and
standard deviation of the untruncated distribution.
"""
function tnvar(a, b, μ, σ)
    α = (a - μ)/σ; β = (b - μ)/σ
    return tnvar(α, β) * σ^2
end
