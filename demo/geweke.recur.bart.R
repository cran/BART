
library(BART)

## simulate recurrent events data set with Exponential proportional intensity
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
    x[k, 4:13] <- runif(10)
    x[k, 14:23] <- rbinom(10, 1, 0.5)

    for(j in 1:K) {
        x[k, 1:3] <- c(j, j-v, N.[k])
        if(j>1) x[k, 4:23] <- x[k-1, 4:23]
        alpha <- 0.0001*exp(sum(b*x[k, 4:23])+sqrt(N.[k]))
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

for(i in 1:N) {
    j <- (i-1)*K+1:K
    if(i==1) plot(1:K, cum[j], type='l',
                  xlab='t', ylab=expression(Lambda(t, x)),
                  sub='Proportional Setting', ylim=c(0, 20))
    else lines(1:K, cum[j], col=i)
}

post <- mc.recur.bart(x.train=x, y.train=y, x.test=x, nskip=1000,
                      keepevery=100, seed=99, mc.cores=C,
                      sparse=TRUE, augment=TRUE)

for(i in 1:N) {
    j <- (i-1)*K+1:K
    if(i==1) plot(1:K, post$cum.test.mean[j], type='l',
                  xlab='t', ylab=expression(Lambda(t, x)),
                  sub='Proportional Setting', ylim=c(0, 20))
    else lines(1:K, post$cum.test.mean[j], col=i)
}

print(cor(cum, post$cum.test.mean)^2)

plot(cum, post$cum.test.mean, pch='.',
     xlim=c(0, 20), ylim=c(0, 20),
     xlab='True', ylab='Estimated', sub='Proportional Setting')
abline(0, 1)

## convergence diagnostics
par(mfrow=c(2, 2))

## select 10 values of x1 for Friedman's partial dependence function
M <- nrow(post$yhat.test)
m <- 10
x1 <- seq(0, 1, length.out=m)
pred <- as.list(1:m)
Fpdf <- as.list(1:m)

y. <- 0
for(i in 1:m) {
    x.test <- x
    x.test[ , 4] <- x1[i]
    if(length(pred[[i]])==1) {
        pred[[i]] <- predict(post, newdata=x.test, mc.cores=C)
        Fpdf[[i]] <- apply(pnorm(pred[[i]]$yhat.test), 1, mean)
    }
    y.[i] <- mean(Fpdf[[i]])
}

plot(x1, y., type='l',
     xlab=expression(x[1]), ylab='partial dependence function')

# select 10 subject X time points to summarize
i <- floor(seq(1, NK, length.out=m))
j <- seq(-0.5, 0.4, length.out=m)
for(h in 1:m) {
    auto.corr <- acf(post$yhat.test[ , i[h]], plot=FALSE)
    if(h==1) {
        max.lag <- max(auto.corr$lag[ , 1, 1])
        plot(1:max.lag+j[h], auto.corr$acf[1+(1:max.lag), 1, 1],
             type='h', xlim=c(0, max.lag+1), ylim=c(-1, 1),
             ylab='auto-correlation', xlab='lag')
    }
    else
        lines(1:max.lag+j[h], auto.corr$acf[1+(1:max.lag), 1, 1],
              type='h', col=h)
}

for(j in 1:m) {
    ##if(j==1) plot(pnorm(post$yhat.test[ , i[j]]), ylim=c(0, 0.16),
    if(j==1) plot(post$yhat.test[ , i[j]], ylim=c(-3, -1),
                  type='l', ylab='f(x)', xlab='m')
    else lines(post$yhat.test[ , i[j]], type='l', col=j)
}

## select 10 subjects uniformly spread out over the data set
h <- seq(1, N*K, floor(N/m)*K)
j <- 1
for(i in h) {
    z <- gewekediag(post$yhat.test[ , (i-1)+1:K])$z
    y <- max(c(4, abs(z)))

    ## plot the z scores vs. time for each patient
    if(i==1) plot(post$times, z, ylim=c(-y, y), type='l',
                  xlab='t', ylab='z')
    else lines(post$times, z, type='l', col=j)
    j <- j+1
}
## add two-sided alpha=0.05 critical value lines
lines(post$times, rep(-1.96, K), type='l', lty=2)
lines(post$times, rep( 1.96, K), type='l', lty=2)

par(mfrow=c(1, 1))
