# TruncatedNormal

[![Build Status](https://travis-ci.org/cossio/TruncatedNormal.jl.svg?branch=master)](https://travis-ci.org/cossio/TruncatedNormal.jl)
[![Coverage Status](https://coveralls.io/repos/cossio/TruncatedNormal.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/cossio/TruncatedNormal.jl?branch=master)
[![codecov.io](http://codecov.io/github/cossio/TruncatedNormal.jl/coverage.svg?branch=master)](http://codecov.io/github/cossio/TruncatedNormal.jl?branch=master)

Install with

    Pkg.clone("https://github.com/cossio/TruncatedNormal.jl")
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

See https://github.com/cossio/TruncatedNormal.jl/blob/master/notes/normal.pdf for mathematical details.