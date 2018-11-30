
library(BART)

N = 1000
P = 5       #number of covariates
ndpost = 1000
nskip = 100
C = 8

set.seed(12)
x.train=matrix(runif(N*P, -2, 2), N, P)
Ey.train = x.train[ , 1]^3
y.train=rbinom(N, 1, pnorm(Ey.train))
table(y.train)

##run BART with C cores in parallel
post = mc.lbart(x.train, y.train, mc.cores=C, keepevery=10,
                seed=99, ndpost=ndpost, nskip=nskip)

M <- 11
x <- seq(-2, 2, length.out=M)
x.test <- as.matrix(expand.grid(list(x5=x, x4=x, x3=x, x2=x, x1=x)))[ , 5:1]
pred <- predict(post, x.test, mc.cores=C)

K <- M^4
prob.mean <- 0
prob.025 <- 0
prob.975 <- 0
for(i in 1:M) {
   j <- (i-1)*K+1:K
   prob <- apply(pred$prob.test[ , j], 1, mean)
   prob.mean[i] <- mean(prob)
   prob.025[i] <- quantile(prob, probs=0.025)
   prob.975[i] <- quantile(prob, probs=0.975)
}

X <- seq(-2, 2, length.out=50)
plot(X, pnorm(X^3), type='l',
     xlab='x', ylab=expression(Phi(f(x))))
lines(x, prob.mean, col='blue')
lines(x, prob.025, col='red')
lines(x, prob.975, col='red')
##dev.copy2pdf(file='cube-lbart.pdf')
