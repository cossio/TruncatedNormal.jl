library(TruncatedNormal)

X <- TruncatedNormal::mvrandn(l = c(724.1281846948646, -0.3248194681994417),
                             u = c(2518.3641949975517, 0.5111929404723518),
                             Sig = matrix(c(1., -0.022020238658485524, -0.022020238658485524, 1.), 2,2),
                             n = 1000000)


tmu <- c(mean(X[1,]), mean(X[2,]))

tSig <- matrix(c(mean(X[1,] * X[1,]), mean(X[1,] * X[2,]), 
                 mean(X[2,] * X[1,]), mean(X[2,] * X[2,])), 2,2)

