
library(BART)

N = 1000
P = 5       #number of covariates
ndpost = 1000
nskip = 100
C = 8

set.seed(12)
x.train=matrix(runif(N*P, -2, 2), N, P)
Ey.train = x.train[ , 1]^3
y.train=rnorm(N, Ey.train)

##run BART with C cores in parallel
post = mc.wbart(x.train, y.train, mc.cores=C,
                seed=99, ndpost=ndpost, nskip=nskip)

gost = mc.gbart(x.train, y.train, mc.cores=C,
                seed=99, ndpost=ndpost, nskip=nskip)

plot(post$yhat.train.mean, gost$yhat.train.mean)
print(cor(post$yhat.train.mean, gost$yhat.train.mean))

set.seed(12)
x.train=matrix(runif(N*P, -2, 2), N, P)
Ey.train = x.train[ , 1]^3
y.train=rbinom(N, 1, pnorm(Ey.train))
table(y.train)

##run BART with C cores in parallel
post = mc.pbart(x.train, y.train, mc.cores=C, keepevery=10,
                seed=99, ndpost=ndpost, nskip=nskip)

gost = mc.gbart(x.train, y.train, mc.cores=C, type='pbart',
                seed=99, ndpost=ndpost, nskip=nskip)

plot(post$prob.train.mean, gost$prob.train.mean,
     xlim=0:1, ylim=0:1)
print(cor(post$prob.train.mean, gost$prob.train.mean))

post = mc.lbart(x.train, y.train, mc.cores=C, keepevery=10,
                seed=99, ndpost=ndpost, nskip=nskip)

gost = mc.gbart(x.train, y.train, mc.cores=C, type='lbart',
                tau.num=3.663562, seed=99, ndpost=ndpost, nskip=nskip)

plot(post$prob.train.mean, gost$prob.train.mean,
     xlim=0:1, ylim=0:1)
print(cor(post$prob.train.mean, gost$prob.train.mean))
