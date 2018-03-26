export tnmean, tnvar

"""
    tnmean(a,b)

Mean of the truncated standard normal distribution

    tnmean(a, b, μ, σ)

Mean of the truncated normal distribution, where μ, σ
are the mean and standard deviation of the untruncated
distribution.
"""
function tnmean end

function tnmean(a, b)
    if a < b
        √(2/π) * _F1(a/√2, b/√2)
    elseif a == b
        a
    else
        throw(ArgumentError("a must be < b"))
    end
end

function tnmean(a, b, μ, σ)
    α = (a - μ)/σ; β = (b - μ)/σ
    return μ + tnmean(α, β) * σ
end


"""
    tnvar(a,b)

Variance of the truncated standard normal distribution

    tnvar(a, b, μ, σ)

Variance of the truncated normal distribution, where μ,σ
are the mean and standard deviation of the untruncated
distribution.
"""
function tnvar end

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

function tnvar(a, b, μ, σ)
    α = (a - μ)/σ; β = (b - μ)/σ
    return tnvar(α, β) * σ^2
end
