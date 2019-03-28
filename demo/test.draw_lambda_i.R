library(BART)

T <- 1
MU <- 0

set.seed(12)

lambda <- draw_lambda_i(1, MU)
rtnorm(1, MU, sqrt(lambda), T)
##rtnorm(MU, T, sqrt(lambda))

set.seed(12)

N <- 10000

lambda <- draw_lambda_i(1, MU)
y <- rtnorm(N, MU, sqrt(lambda), T)
##y <- rtnorm(MU, T, sqrt(lambda))

for(i in 2:N) {
    lambda[i] <- draw_lambda_i(lambda[i-1], MU)
    ##y[i] <- rtnorm(MU, T, sqrt(lambda[i]))
}

x <- seq(T, T+2, length.out=1000)

plot(x, dlogis(x, MU, 1)/plogis(T, MU, 1, lower.tail=FALSE),
     lty=2, type='l',
     ylab=expression(Logistic(MU, 1)))
lines(density(y))
abline(v=T)

##dev.copy2pdf(file='test.draw_lambda_i.pdf')
