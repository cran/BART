
library(BART)

## simulate recurrent events data set with Exponential nonproportional intensity
N <- 250
K <- 60
NK <- N*K
C <- 8

set.seed(-1)

x <- matrix(nrow=NK, ncol=23)
dimnames(x)[[2]] <- c('t', 'v', 'N', paste0('x', 1:20))

b <- c(1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       1.5, 0, 0, 0, 0, 0, 0, 0, 0, 0)

N. <- double(NK)
y <- integer(NK)
cum <- double(NK)
k <- 1

for(i in 1:N) {
    v <- 0
    x[k,  4:13] <- runif(10)
    x[k, 14:23] <- rbinom(10, 1, 0.5)

    for(j in 1:K) {
        x[k, 1:3] <- c(j, j-v, N.[k])
        if(j>1) x[k, 4:23] <- x[k-1, 4:23]
        alpha <- 0.0001*exp(sum(b*x[k, 4:23])*2*(N.[k]+1)/sqrt(j))
        cum[k] <- pexp(30, alpha)
        y[k] <- rbinom(1, 1, cum[k])

        if(y[k]==1) v <- j

        if(j>1) cum[k] <- cum[k-1]+cum[k]
        if(j<K) N.[k+1] <- N.[k]+y[k]

        k <- k+1
    }
}

table(x[K*(1:N), 3])
table(x[K*(1:N), 3])/N

post <- mc.recur.bart(x.train=x, y.train=y, x.test=x,
                      nskip=1000, keepevery=100,
                      sparse=TRUE, mc.cores=C, seed=99)

print(cor(cum, post$cum.test.mean)^2)

par(mfrow=c(1, 2))

for(i in 1:N) {
    j <- (i-1)*K+1:K
    if(i==1) plot(1:K, cum[j], type='l',
                  xlab='t', ylab=expression(Lambda(t, x)),
                  sub='Known Values', ylim=c(0, 20))
    else lines(1:K, cum[j], col=i)
}

for(i in 1:N) {
    j <- (i-1)*K+1:K
    if(i==1) plot(1:K, post$cum.test.mean[j], type='l',
                  xlab='t', ylab=expression(Lambda(t, x)),
                  sub='Estimates', ylim=c(0, 20))
    else lines(1:K, post$cum.test.mean[j], col=i)
}

par(mfrow=c(1, 1))

dev.copy2pdf(file='np-recur-bart.pdf')
