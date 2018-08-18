
library(BART)

set.seed(12)
N=50
P=3

x.train=matrix(runif(N*P, -1, 1), nrow=N, ncol=P)
y=x.train[ , 1]^3
x.miss=matrix(1*(runif(N*P)<0.05), nrow=N, ncol=P)
x.train=x.train*(1-x.miss)
x.train[x.train==0]=NA

post=gbart(x.train, y, x.train)

summary(post$yhat.train.mean)
summary(post$yhat.test.mean)

plot(post$yhat.train.mean, post$yhat.test.mean)

