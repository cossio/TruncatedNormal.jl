library(TruncatedNormal)

X = TruncatedNormal::mvrandn(l = c(-1,1), u = c(-1,1), 
                             Sig =  matrix(c(1, 0.9, 0.9, 1), 2,2),
                             n = 100)
