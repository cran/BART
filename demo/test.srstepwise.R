
library(BART)

set.seed(12)
N=500
P=501
X=matrix(runif(N*P, -1, 1), nrow=N, ncol=P)
dimnames(X)[[2]]=paste0('x', 1:P)
y=rnorm(N, (X[ , 1]^3)+(X[ , 2]^3)+(X[ , 3]^3)+(X[ , 4]^3)+(X[ , 5]^3))
T=exp(y)
C=rexp(N, 0.65)
delta=(T<C)*1
table(delta)/N
times=T*delta+C*(1-delta)

check=srstepwise(X, times, delta)
print(check)
