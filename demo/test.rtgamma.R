library(BART)

N <- 1
A <- 3
SHAPE <- 5
RATE <- 0.5

set.seed(12)

rtgamma(N, SHAPE, RATE, A)

set.seed(12)

rtgamma(N, SHAPE, RATE, A)

set.seed(12)

N <- 10000

y <- 0

y <- rtgamma(N, SHAPE, RATE, A)
##for(i in 1:N) y[i] <- rtgamma(SHAPE, RATE, A)

x <- seq(A, 4*A, length.out=1000)
plot(x, dgamma(x, SHAPE, RATE)/pgamma(A, SHAPE, RATE, lower.tail=FALSE),
     lty=2, type='l', ylim=c(0, 1),
     ylab=expression(Gam(SHAPE, RATE)))
lines(density(y, from=A), col='red')
abline(v=A)
dev.copy2pdf(file='test.rtgamma.pdf')
