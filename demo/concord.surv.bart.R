
library(BART)

N <- 1000
K <- 25
set.seed(12)
x <- rbinom(N, 1, 0.5)
l1 <- 0.5
l2 <- 2
T <- rexp(N, l1*(1-x)+x*l2)
C <- rexp(N, 0.25)
delta <- (T < C)
times <- delta*T+(1-delta)*C
table(delta)/N

post <- mc.surv.bart(x.train=x, times=times, delta=delta,
                     x.test=matrix(0:1, nrow=2, ncol=1),
                     K=25, mc.cores=8, seed=99)

c.true <- l1/(l1+l2)
c.true
c.est <- (1-post$prob.test[ , 1])-(1-post$prob.test[ , 26])
for(j in 2:K)
    c.est <- c.est+((1-post$prob.test[ , j])-(1-post$prob.test[ , j+K]))*
        post$surv.test[ , j-1]*post$surv.test[ , j-1+K]
c.est <- 0.5*(1-c.est)
mean(c.est)
quantile(c.est, probs=c(0.025, 0.975))
mean(c.est-c.true)
quantile(c.est-c.true, probs=c(0.025, 0.975))
