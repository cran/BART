
library(BART)

k = 50       #number of covariates
C = 8
thin = 10
ndpost = 1000
nskip = 100

##pdf(file='geweke-pbart1.pdf')
par(mfrow=c(2, 2))

for(n in c(100, 1000, 10000)) {
    set.seed(12)
    x.train=matrix(runif(n*k, min=-1, max=1), n, k)
    Ey.train = pnorm(x.train[ , 1]^3)
    y.train=rbinom(n, 1, Ey.train)
    table(y.train)

    ##run BART with C cores in parallel
    mc.train = mc.pbart(x.train, y.train, mc.cores=C, keepevery=thin,
                        seed=99, ndpost=ndpost, nskip=nskip)

    x <- seq(-1, 1, length.out=200)

    plot(x, x^3, ylab='f(x)', type='l',
         sub=paste0('N:', n, ', k:', k, ', thin:', thin))
    points(x.train[ , 1], mc.train$yhat.train.mean, pch='.', cex=2)

    i <- floor(seq(1, n, length.out=10))

    auto.corr <- acf(mc.train$yhat.train[ , i], plot=FALSE)
    max.lag <- max(auto.corr$lag[ , 1, 1])

    j <- seq(-0.5, 0.4, length.out=10)
    for(h in 1:10) {
        if(h==1)
            plot(1:max.lag+j[h], auto.corr$acf[1+(1:max.lag), h, h],
                 type='h', xlim=c(0, max.lag+1), ylim=c(-1, 1),
                 ylab='acf', xlab='lag')
        else
            lines(1:max.lag+j[h], auto.corr$acf[1+(1:max.lag), h, h],
                 type='h', col=h)
    }

    for(j in 1:10) {
        if(j==1)
            plot(pnorm(mc.train$yhat.train[ , i[j]]),
                 type='l', ylim=c(0, 1),
                 sub=paste0('N:', n, ', k:', k, ', thin:', thin),
                 ylab=expression(Phi(f(x))), xlab='m')
        else
            lines(pnorm(mc.train$yhat.train[ , i[j]]),
                 type='l', col=j)
    }

    geweke <- gewekediag(mc.train$yhat.train)

    j <- -10^(log10(n)-1)
    plot(geweke$z, pch='.', cex=2, ylab='z', xlab='i',
         sub=paste0('N:', n, ', k:', k, ', thin:', thin),
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
}
##dev.off()

par(mfrow=c(1, 1))
