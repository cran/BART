library(BART)

T <- 8
MU <- 5
SD <- 0.5

set.seed(12)

rtnorm(MU, T, SD)

set.seed(12)

rtnorm(MU, T, SD)

set.seed(12)

N <- 10000

y <- 0

for(i in 1:N) y[i] <- rtnorm(MU, T, SD)

x <- seq(T, T+2*SD, length.out=1000)

plot(x, dnorm(x, MU, SD)/pnorm(T, MU, SD, lower.tail=FALSE),
     lty=2, type='l',
     ylab=expression(N(MU, SD^2)))
lines(density(y))
abline(v=T)

dev.copy2pdf(file='test.rtnorm.pdf')
