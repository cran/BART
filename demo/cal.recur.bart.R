
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

M <- 5
H <- NK/5

fits <- as.list(1:M)
x.test <- as.list(1:M)
x.train <- as.list(1:M)
y.train <- as.list(1:M)
xinfo <- bartModelMatrix(x, numcut=100)$xinfo

for(m in 1:M) {
    h <- (m-1)*H+1:H
    x.test[[m]] <- x[h, ]
    x.train[[m]] <- x[-h, ]
    y.train[[m]] <- y[-h]
    fits[[m]] <- mc.recur.bart(x.train=x.train[[m]], y.train=y.train[[m]],
                               x.test=x.test[[m]], nskip=1000, keepevery=100,
                               sparse=TRUE, augment=TRUE,
                               keeptrainfits=TRUE, xinfo=xinfo, mc.cores=8,
                               seed=m)
}

for(m in 1:M) {
    h <- (m-1)*H+1:H
    if(m==1) {
        cum.test.mean <- fits[[1]]$cum.test.mean
        cum.test <- fits[[1]]$cum.test
        cum.train <- cbind(fits[[M]]$cum.train[ , h], fits[[1]]$cum.train)
    }
    else {
        k <- 1:(m*H)
        cum.test.mean <- c(cum.test.mean, fits[[m]]$cum.test.mean)
        cum.test <- cbind(cum.test, fits[[m]]$cum.test)
        if(m<M) cum.train <- rbind(cum.train, cbind(fits[[m]]$cum.train[ , k],
                                                    fits[[M]]$cum.train[ , h],
                                                    fits[[m]]$cum.train[ , -k]))
    }
}

## equal tailed
lower <- 0
upper <- 1
tol <- 0.001
cover <- 1
value <- 0.95
while(!((value-tol)<=cover & cover<=(value+tol))) {
print(c(cover=cover, lower=lower, upper=upper))
cum.train.lower <- apply(cum.train, 2, quantile, probs=lower)
cum.train.upper <- apply(cum.train, 2, quantile, probs=upper)
cover <- mean(cum.train.lower<=cum.test.mean & cum.test.mean<=cum.train.upper)
if(cover<value) {
    upper <- upper+0.01
    if(upper>1) break
}
else upper <- upper/1.01
lower <- 1-upper
}
print(c(cover=cover, lower=lower, upper=upper))


