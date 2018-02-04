
library(BART)

## estimate concordance probability: P(t1<t2)
N <- 2000
r1 <- 0.1
r2 <- 0.3
## true concordance
true.p <- r1/(r1+r2)

set.seed(12)
x <- rbinom(N, 1, 0.5)
t <- ceiling(rexp(N, r1*x+r2*(1-x)))
c <- ceiling(rexp(N, 0.035))
delta <- (t<c)
table(delta)/N
t <- delta*t+(1-delta)*c

post <- mc.surv.bart(x.train=cbind(x), x.test=rbind(0, 1),
                     times=t, delta=delta, mc.cores=8, seed=99) 

K <- post$K
q <- 1-pnorm(post$yhat.test)
for(j in 1:K) {
    if(j==1) P <- q[ , K+1]-q[ , 1]
    else P <- P+(q[ , K+j]-q[ , j])*post$surv.test[ , K+j-1]*post$surv.test[ , j-1]
}
C <- 0.5*(1-P)
## estimate concordance
summary(C)
