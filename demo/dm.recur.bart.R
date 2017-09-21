
library(BART)

## load 20 percent random sample
data(xdm20.train)
data(xdm20.test)
data(ydm20.train)

set.seed(99)
post <- recur.bart(x.train=xdm20.train, y.train=ydm20.train)
## larger data sets can take some time so, if parallel processing
## is available, submit this statement instead
## post <- mc.recur.bart(x.train=xdm20.train, y.train=ydm20.train,
##                      mc.cores=8, seed=99)

require(rpart)
require(rpart.plot)

dss <- rpart(post$yhat.train.mean~xdm20.train)

rpart.plot(dss)
## for the 20 percent sample, notice that the top splits
## involve cci_pvd and n
## for the full data set, notice that all splits
## involve ca, cci_pud, cci_pvd, ins270 and n
## (except one at the bottom involving a small group)

## compare patients treated with insulin (ins270=1) vs
## not treated with insulin (ins270=0)
N.train <- 50
N.test <- 50
K <- post$K ## 798 unique time points

## only testing set, i.e., remove training set
xdm20.test. <- xdm20.test[N.train*K+(1:(N.test*K)), ]
xdm20.test. <- rbind(xdm20.test., xdm20.test.)
xdm20.test.[ , 'ins270'] <- rep(0:1, each=N.test*K)

## multiple threads will be utilized if available
pred <- predict(post, xdm20.test., mc.cores=8)

## create Friedman's partial dependence function for the
## intensity/hazard by time and ins270
NK.test <- N.test*K
M <- nrow(pred$haz.test) ## number of MCMC samples, typically 1000

RI <- matrix(0, M, K)

for(i in 1:N.test)
    RI <- RI+(pred$haz.test[ , (N.test+i-1)*K+1:K]/
              pred$haz.test[ , (i-1)*K+1:K])/N.test

RI.lo <- apply(RI, 2, quantile, probs=0.025)
RI.mu <- apply(RI, 2, mean)
RI.hi <- apply(RI, 2, quantile, probs=0.975)

plot(post$times, RI.hi, type='l', lty=2, log='y',
     ylim=c(min(RI.lo, 1/RI.hi), max(1/RI.lo, RI.hi)),
     xlab='t', ylab='RI(t, x)',
     sub='insulin(ins270=1) vs. no insulin(ins270=0)',
     main='Relative intensity of hospital admissions for diabetics')
lines(post$times, RI.mu)
lines(post$times, RI.lo, lty=2)
lines(post$times, rep(1, K), col='darkgray')

## RI for insulin therapy seems fairly constant with time
mean(RI.mu)
