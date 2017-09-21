
library(BART)

n = 1000
k = 5       #number of covariates
ndpost = 1000
nskip = 100
C = 8

set.seed(12)
x.train=matrix(runif(n*k, -2, 2), n, k)
Ey.train = x.train[ , 1]^3
y.train=rbinom(n, 1, pnorm(Ey.train))
table(y.train)

##run BART with C cores in parallel
mc.train = mc.pbart(x.train, y.train, mc.cores=C, keepevery=10,
                    seed=99, ndpost=ndpost, nskip=nskip)

yhat.train.025 <- apply(mc.train$yhat.train, 2, quantile, probs=0.025)
yhat.train.975 <- apply(mc.train$yhat.train, 2, quantile, probs=0.975)

x <- seq(-2, 2, length.out=200)

plot(x, x^3, type='l', ylab='f(x)')
points(x.train[ , 1], mc.train$yhat.train.mean, pch='.', cex=2)
