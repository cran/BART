
library(BART)

N = 500
P = 1       #number of covariates
C = 8
SD.y=10
M = 201

set.seed(12)
x.train=matrix(runif(N*P, -2, 2), N, P)
x.test=matrix(seq(-2, 2, length.out=M), M, P)
Ey.train = x.train[ , 1]^3
y.train=rnorm(N, Ey.train, SD.y)

post <- as.list(1:4)

post[[1]] = mc.gbart(x.train, y.train, x.test, mc.cores=C, seed=99)

info <- c(0, N, 4*N)

for(i in 2:3) {
    M=info[i]
    x.info=matrix(runif(M*P, -2, 2), M, P)
    Ey.info = x.info[ , 1]^3
    y.info=c(rnorm(M, Ey.info, SD.y), y.train)
    x.info=rbind(x.info, x.train)
    post[[i]] = mc.gbart(x.info, y.info, x.test, mc.cores=C, seed=99)
}

## plot(x.test[ , 1], x.test[ , 1]^3, type='l', xlab='x', ylab='y')
## legend(-2, 8, lty=rep(1:0, 4), legend=c('Truth', ' ', 'No prior', ' ',
##                               'Equiv. N', ' ', 'Equiv. 4N', ' '),
##        col=rep(1:4, each=2), cex=0.8, bty='n')
## legend(-0.75, 8, lty=rep(0, 8),
##        legend=c(expression(italic(R)^2), ' ', '0.828', ' ',
##                 '0.858', ' ', '0.881', ' '),
##        col=rep(1:4, each=2), cex=0.8, bty='n')
## legend(0, 8, lty=rep(0, 8), legend=c('MSE', ' ', '1.93', ' ',
##                                      '1.73', ' ', '1.24', ' '),
##        col=rep(1:4, each=2), cex=0.8, bty='n')

## for(i in 1:3) {
##     print(cor(x.test[ , 1]^3, post[[i]]$yhat.test.mean)^2)
##     lines(x.test[ , 1], post[[i]]$yhat.test.mean, col=i+1)
## }
## dev.copy2pdf(file='inform-alt.pdf', height=4, width=4)

for(i in 1:3) 
    print(mean(((x.test[ , 1]^3)-post[[i]]$yhat.test.mean)^2))
