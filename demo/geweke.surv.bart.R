
library(BART)

f <- function(x) ## adpated Friedman's five dimensional test function
    3+sin(pi*x[ , 1]*x[ , 2])-2*(x[ , 3]-0.5)^2+x[ , 4]-0.5*x[ , 5]

N <- 100
A <- 1/155
P <- 20 #number of covariates

set.seed(12)

x.train <- matrix(runif(N*P), nrow=N, ncol=P)
x.test <- matrix(runif(N*P), nrow=N, ncol=P)

T <- ceiling(rweibull(N, shape=2, scale=exp(f(x.train))))
C <- ceiling(rexp(N, A))
delta <- 1*(T<C)
times <- delta*T+(1-delta)*C
table(delta)/N

C = 8

##run BART with C cores in parallel
post = mc.surv.bart(x.train, times=times, delta=delta, mc.cores=C, seed=99,
                    keepevery=50)

x <- x.train

x4 <- seq(0, 1, length.out=10)

for(i in 1:10) {
    x[ , 4] <- x4[i]

    if(i==1) x.test <- x
    else x.test <- rbind(x.test, x)
}

pre = surv.pre.bart(times=times, delta=delta,
                    x.train=x.train, x.test=x.test)

##run predict with C cores in parallel
pred <- predict(post, newdata=pre$tx.test, mc.cores=C)

K <- pred$K

##create Friedman's partial dependence function for x4
surv <- list(1:10)

for(i in 1:10) {
    surv[[i]] <- matrix(nrow=1000, ncol=K)

    for(j in 1:K) {
        h <- (i-1)*N*K+seq(j, N*K, by=K)
        surv[[i]][ , j] <- apply(pred$surv.test[ , h], 1, mean)
    }

    surv[[i]] <- apply(surv[[i]], 2, mean)
}

for(i in 1:10) {
    if(i==1) plot(c(0, pre$times), c(1, surv[[i]]), type='s',
                  xlim=c(0, 50), ylim=0:1, xlab='t', ylab='S(t, x)')
    else lines(c(0, pre$times), c(1, surv[[i]]), type='s', col=i)
    j <- min(which(surv[[i]]<0.5))
    text(pre$times[j], 0.5, paste(round(x4[i], digits=1)), col=i, pos=2)
}

## acf plots for 10 subjects
k <- floor(seq(1, N, length.out=10))
j. <- seq(-0.5, 0.4, length.out=10)

for(j in 1:K) {
    for(i in 1:10) {
        h <- (k[i]-1)*K+j

        auto.corr <- acf(pred$yhat.test[ , h], plot=FALSE)
        max.lag <- max(auto.corr$lag[ , 1, 1])

            if(i==1)
                plot(1:max.lag+j.[i], auto.corr$acf[1+(1:max.lag), 1, 1],
                     type='h', xlim=c(0, max.lag+1), ylim=c(-1, 1),
                     sub=paste0('t=', pre$times[j]), ylab='acf', xlab='lag')
            else
                lines(1:max.lag+j.[i], auto.corr$acf[1+(1:max.lag), 1, 1],
                      type='h', col=i)
        }

    Sys.sleep(1)
}

## trace plots for 10 subjects
k <- floor(seq(1, N, length.out=10))

for(j in 1:K) {
    for(i in 1:10) {
                   
    h <- (k[i]-1)*K+j

    if(i==1)
        plot(pred$yhat.test[ , h], type='l',
             ylim=c(-4, 0), sub=paste0('t=', pre$times[j]),
             ylab=expression(Phi(f(x))), xlab='m')
    else
        lines(pred$yhat.test[ , h], type='l', col=i)
    }
    Sys.sleep(1)
}

## Geweke plot for 10 subjects
k <- floor(seq(1, N, length.out=10))

geweke <- list(1:10)

for(i in 1:10) {
    h <- (k[i]-1)*K+1:K
    
    geweke[[i]] <- gewekediag(pred$yhat.test[ , h])
}

max.t <- max(pre$times)
min.t <- -max.t/10

for(i in 1:10) {
    if(i==1) {
        plot(pre$times, geweke[[i]]$z, type='l',
             ylab='z', xlab='t', ylim=c(-5, 5), xlim=c(min.t, max.t))
        lines(pre$times, rep(-1.96, K), type='l', col=6)
        lines(pre$times, rep(+1.96, K), type='l', col=6)
        lines(pre$times, rep(-2.576, K), type='l', col=5)
        lines(pre$times, rep(+2.576, K), type='l', col=5)
        lines(pre$times, rep(-3.291, K), type='l', col=4)
        lines(pre$times, rep(+3.291, K), type='l', col=4)
        lines(pre$times, rep(-3.891, K), type='l', col=3)
        lines(pre$times, rep(+3.891, K), type='l', col=3)
        lines(pre$times, rep(-4.417, K), type='l', col=2)
        lines(pre$times, rep(+4.417, K), type='l', col=2)
        text(c(0, 0), c(-1.96, 1.96), pos=2, cex=0.6, labels='0.95')
        text(c(0, 0), c(-2.576, 2.576), pos=2, cex=0.6, labels='0.99')
        text(c(0, 0), c(-3.291, 3.291), pos=2, cex=0.6, labels='0.999')
        text(c(0, 0), c(-3.891, 3.891), pos=2, cex=0.6, labels='0.9999')
        text(c(0, 0), c(-4.417, 4.417), pos=2, cex=0.6, labels='0.99999')
    }
    else lines(pre$times, geweke[[i]]$z, type='l')
}
