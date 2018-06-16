library(TruncatedNormal)

u = c(-0.025874, 0.0414868)
Sig = matrix(c(4.23262e-7, 2.95285e-6, 2.95285e-6, 2.06834e-5), 2, 2)
lb = c(0,0)
ub = c(1000, 0.0970138)

X <- TruncatedNormal::mvrandn(l = lb - u, u = ub - u, Sig = Sig, n = 1000000)
X <- X + u

tmu <- c(mean(X[1,]), mean(X[2,]))

tSig <- matrix(c(mean(X[1,] * X[1,]), mean(X[1,] * X[2,]), 
                 mean(X[2,] * X[1,]), mean(X[2,] * X[2,])), 2,2)

