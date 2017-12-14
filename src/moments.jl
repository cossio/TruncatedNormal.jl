export tnmean, tnvar

"""
    tnmean(a,b)
Mean of the truncated standard normal distribution

    tnmean(a, b, μ, σ)
Mean of the truncated normal distribution, where μ,σ
are the mean and standard deviation of the untruncated
distribution.
"""
function tnmean end

function tnmean(a, b)
    if -Inf < a ≤ b < Inf
        √(2/π) * _F1(a/√2, b/√2)
    else
        throw(ArgumentError("tnmean called with invalid parameters a=$a, b=$b"))
    end
end

function tnmean(a, b, μ, σ)
    if -Inf < a ≤ b < Inf && -Inf < μ < Inf && 0 < σ < Inf
        α = (a - μ)/σ; β = (b - μ)/σ
        μ + tnmean(α, β) * σ
    else
        throw(ArgumentError("tnmean called with invalid parameters a=$a, b=$b, μ=$μ, σ=$σ"))
    end
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
    if -Inf < a ≤ b < Inf
        1 + 2/√π * _F2(a/√2, b/√2) - 2/π * _F1(a/√2, b/√2)^2
    else
        throw(ArgumentError("tnvar called with invalid parameters a=$a, b=$b"))
    end
end

function tnvar(a, b, μ, σ)
    if -Inf < a ≤ b < Inf && -Inf < μ < Inf && 0 < σ < Inf
        α = (a - μ)/σ; β = (b - μ)/σ
        tnvar(α, β) * σ^2
    else
        throw(ArgumentError("tnvar called with invalid parameters a=$a, b=$b, μ=$μ, σ=$σ"))
    end
end
