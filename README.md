# TruncatedNormal

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://cossio.github.io/TruncatedNormal.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://cossio.github.io/TruncatedNormal.jl/dev)
![](https://github.com/cossio/TruncatedNormal.jl/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/cossio/TruncatedNormal.jl/branch/master/graph/badge.svg?token=c1Qv0pcqn5)](https://codecov.io/gh/cossio/TruncatedNormal.jl)

Install with

    Pkg.add("https://github.com/cossio/TruncatedNormal.jl")
    using TruncatedNormal

Mean of the truncated standard normal distribution:

    tnmean(a,b)

Mean of the truncated normal distribution, where μ,σ
are the mean and standard deviation of the untruncated
distribution:

    tnmean(a, b, μ, σ)

Variance of the truncated standard normal distribution:

    tnvar(a,b)

Variance of the truncated normal distribution, where μ,σ
are the mean and standard deviation of the untruncated
distribution:

    tnvar(a, b, μ, σ)

It works even if the truncation interval is far from the mode of the distribution. See https://github.com/cossio/TruncatedNormal.jl/blob/master/notes/normal.pdf for mathematical details.
