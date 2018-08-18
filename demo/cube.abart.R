
library(BART)

N = 1000
NP = 200
P = 5       #number of covariates
M = 8
ndpost = 1000
set.seed(12)
x.train=matrix(runif(N*P, -2, 2), N, P)
mu = x.train[ , 1]^3
print(quantile(mu, probs=c(0.1, 0.3, 0.5, 0.7, 0.9)))
x.test=matrix(runif(NP*P, -2, 2), NP, P)
x1=c(-1, 0, 1)
x.test=cbind(c(rep(x1[1], NP), rep(x1[2], NP), rep(x1[3], NP)),
             rbind(x.test, x.test, x.test)[ , -1])

y=rnorm(N, mu)
offset=mean(y)
T=exp(y)
C=rexp(N, 0.05)
delta=(T<C)*1
table(delta)/N
times=(T*delta+C*(1-delta))

post = mc.abart(x.train, times, delta, x.test,
                 mc.cores=M, seed=99, ndpost=ndpost)
## post2 = mc.abart(x.train, times, delta, x.test, offset=offset,
##                  mc.cores=M, seed=99, ndpost=ndpost)

Z=8

plot(mu, post$yhat.train.mean, asp=1,
     xlim=c(-Z, Z), ylim=c(-Z, Z))
abline(a=0, b=1)

## plot(mu, post2$yhat.train.mean, asp=1,
##      xlim=c(-Z, Z), ylim=c(-Z, Z))
## abline(a=0, b=1)

## plot(post$yhat.train.mean, post2$yhat.train.mean, asp=1,
##      xlim=c(-Z, Z), ylim=c(-Z, Z))
## abline(a=0, b=1)

K <- post$K
par(mfrow=c(3, 1))
for(i in 1:length(x1)) {
        plot(c(0, post$times),
             c(1, pnorm(log(post$times), mean=x1[i]^3,
                   lower.tail=FALSE)),
             type='l', ylim=0:1, xlab='t', ylab='S(t, x)')
        lines(c(0, post$times),
              c(1, post$surv.test.mean[(i-1)*K+1:K]),
              col=2, type='s')
        post$surv.test.025 <- matrix(nrow=ndpost, ncol=K)
        post$surv.test.975 <- matrix(nrow=ndpost, ncol=K)
        for(j in 1:K) {
            post$surv.test.025[ , j] <-
                apply(post$surv.test[ , (i-1)*K*NP+seq(j, K*NP, K)],
                      1, mean)
            post$surv.test.975[ , j] <-
                apply(post$surv.test[ , (i-1)*K*NP+seq(j, K*NP, K)],
                      1, mean)
        }
        post$surv.test.025 <- apply(post$surv.test.025, 2,
                                    quantile, probs=0.025)
        post$surv.test.975 <- apply(post$surv.test.975, 2,
                                    quantile, probs=0.975)
        lines(c(0, post$times), c(1, post$surv.test.025),
              col=2, type='s', lty=2)
        lines(c(0, post$times), c(1, post$surv.test.975),
              col=2, type='s', lty=2)
}
par(mfrow=c(1, 1))
