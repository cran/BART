
library(BART)

N = 1000
P = 5       #number of covariates
C = 8

set.seed(12)
x.train=matrix(runif(N*P, -2, 2), N, P)
Ey.train = x.train[ , 1]^3
y.train=rnorm(N, Ey.train, 3)

##run BART with C cores in parallel
post.est = mc.gbart(x.train, y.train, mc.cores=C, seed=99,
                    sigest=3)

post.fix = mc.gbart(x.train, y.train, mc.cores=C, seed=99,
                    lambda=0, sigest=3)

plot(post.est$yhat.train.mean, post.fix$yhat.train.mean)
print(cor(post.est$yhat.train.mean, post.fix$yhat.train.mean))

