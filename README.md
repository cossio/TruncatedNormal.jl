# TruncatedNormal

**Note:** Merged in to Distributions.jl, see JuliaStats/Distributions.jl#1058, JuliaStats/Distributions.jl#691.

## Installation

Install with

```julia
using Pkg
Pkg.add(url="https://github.com/cossio/TruncatedNormal.jl")
```

## Usage

This package does not export any symbols. The following functions are defined.

Mean of the truncated standard normal distribution:

```julia
TruncatedNormal.tnmean(a,b)
```

Mean of the truncated normal distribution, where μ,σ
are the mean and standard deviation of the untruncated
distribution:

```julia
TruncatedNormal.tnmean(a, b, μ, σ)
```

Variance of the truncated standard normal distribution:

```julia
TruncatedNormal.tnvar(a,b)
```

Variance of the truncated normal distribution, where μ,σ
are the mean and standard deviation of the untruncated
distribution:

```julia
TruncatedNormal.tnvar(a, b, μ, σ)
```

It works even if the truncation interval is far from the mode of the distribution.
See https://github.com/cossio/TruncatedNormal.jl/blob/master/notes/normal.pdf for mathematical details.
