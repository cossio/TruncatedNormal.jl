tnmean(a, b) = √(2/π) * F1(a/√2, b/√2)
tnmean(a, b, μ, σ) = μ + tnmean(a,b) * σ

tnvar(a,b) = 1 + 2/√π * F2(a/√2, b/√2) - 2/π * F1(a/√2, b/√2)^2
tnvar(a,b,μ,σ) = tnvar(a,b) * σ^2
