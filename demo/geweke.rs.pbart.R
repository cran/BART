
library(BART)

##simulate from Friedman's five-dimensional test function
##Friedman JH. Multivariate adaptive regression splines
##(with discussion and a rejoinder by the author).
##Annals of Statistics 1991; 19:1-67.

f = function(x) #only the first 5 matter
    sin(pi*x[ , 1]*x[ , 2]) + 2*(x[ , 3]-.5)^2+x[ , 4]+0.5*x[ , 5]-1.5

sigma = 1.0  #y = f(x) + sigma*z where z~N(0, 1)
k = 50       #number of covariates
thin = 25
ndpost = 2500
nskip = 100
C = 10
m = 10
n = 10000

set.seed(12)
x.train=matrix(runif(n*k), n, k)
Ey.train = f(x.train)
y.train=(Ey.train+sigma*rnorm(n)>0)*1
table(y.train)/n

x <- x.train
x4 <- seq(0, 1, length.out=m)

for(i in 1:m) {
    x[ , 4] <- x4[i]

    if(i==1) x.test <- x
    else x.test <- rbind(x.test, x)
}

post = rs.pbart(x.train, y.train, x.test=x.test,
                C=C, mc.cores=8, keepevery=thin,
                seed=99, ndpost=ndpost, nskip=nskip)
str(post)

par(mfrow=c(2, 2))

M <- nrow(post$yhat.test)
pred <- matrix(nrow=M, ncol=10)

for(i in 1:m) {
    h <- (i-1)*n+1:n
    pred[ , i] <- apply(pnorm(post$yhat.test[ , h]), 1, mean)
}

pred <- apply(pred, 2, mean)

plot(x4, qnorm(pred), xlab=expression(x[4]),
     ylab='partial dependence function', type='l')

i <- floor(seq(1, n, length.out=10))
j <- seq(-0.5, 0.4, length.out=10)
for(h in 1:10) {
    auto.corr <- acf(post$yhat.shard[ , i[h]], plot=FALSE)
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

for(j in 1:10) {
    if(j==1)
        plot(pnorm(post$yhat.shard[ , i[j]]),
             type='l', ylim=c(0, 1),
             sub=paste0('N:', n, ', k:', k),
             ylab=expression(Phi(f(x))), xlab='m')
    else
        lines(pnorm(post$yhat.shard[ , i[j]]),
              type='l', col=j)
}

geweke <- gewekediag(post$yhat.shard)

j <- -10^(log10(n)-1)
plot(geweke$z, pch='.', cex=2, ylab='z', xlab='i',
     sub=paste0('N:', n, ', k:', k),
     xlim=c(j, n), ylim=c(-5, 5))
lines(1:n, rep(-1.96, n), type='l', col=6)
lines(1:n, rep(+1.96, n), type='l', col=6)
lines(1:n, rep(-2.576, n), type='l', col=5)
lines(1:n, rep(+2.576, n), type='l', col=5)
lines(1:n, rep(-3.291, n), type='l', col=4)
lines(1:n, rep(+3.291, n), type='l', col=4)
lines(1:n, rep(-3.891, n), type='l', col=3)
lines(1:n, rep(+3.891, n), type='l', col=3)
lines(1:n, rep(-4.417, n), type='l', col=2)
lines(1:n, rep(+4.417, n), type='l', col=2)
text(c(1, 1), c(-1.96, 1.96), pos=2, cex=0.6, labels='0.95')
text(c(1, 1), c(-2.576, 2.576), pos=2, cex=0.6, labels='0.99')
text(c(1, 1), c(-3.291, 3.291), pos=2, cex=0.6, labels='0.999')
text(c(1, 1), c(-3.891, 3.891), pos=2, cex=0.6, labels='0.9999')
text(c(1, 1), c(-4.417, 4.417), pos=2, cex=0.6, labels='0.99999')

par(mfrow=c(1, 1))

##dev.copy2pdf(file='geweke.rs.pbart.pdf')
