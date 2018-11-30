
library(BART)

##simulate from Friedman's five-dimensional test function
##Friedman JH. Multivariate adaptive regression splines
##(with discussion and a rejoinder by the author).
##Annals of Statistics 1991; 19:1-67.

f = function(x) #only the first 5 matter
    sin(pi*x[ , 1]*x[ , 2]) + 2*(x[ , 3]-.5)^2+x[ , 4]+0.5*x[ , 5]-1.5

sigma = 1.0  #y = f(x) + sigma*z where z~N(0, 1)
k = 100       #number of covariates
ndpost = 1000
nskip = 100
C = 8

par(mfrow=c(3, 2))

post <- as.list(1:6)

for(i in 1:3) {
    n <- 10^(1+i)
    set.seed(12)
    x.train=matrix(runif(n*k), n, k)
    Ey.train = f(x.train)
    y.train=Ey.train+sigma*rnorm(n)

    for(j in c(TRUE, FALSE)) {
        h <- (i-1)*2+j+1
        post[[h]] = mc.wbart(x.train, y.train, mc.cores=C, sparse=TRUE,
                             augment=j, seed=99, ndpost=ndpost, nskip=nskip)

        plot(post[[h]]$varprob.mean, col=c(rep(2, 5), rep(1, k-5)),
             main=paste0('N:', n, ', P:', k, ', Assumption:', c(2.2, 2.1)[j+1]),
             ##sub=expression(-1.5+sin(pi*x[1]*x[2]) + 2*(x[3]-.5)^2+x[4]+0.5*x[5]),
             ylab='Selection Probability', ylim=0:1)
        lines(c(0, 100), c(1/k, 1/k))
    }
}

par(mfrow=c(1, 1))

##dev.copy2pdf(file='sparse.wbart.pdf')

## check1=pwbart(x.train, post[[6]]$treedraws, post[[6]]$mu)
## plot(apply(check1, 2, mean), post[[6]]$yhat.train.mean)

## check2=pwbart(x.train, post[[6]]$treedraws, post[[6]]$mu, dodraws=FALSE)
## plot(check2, post[[6]]$yhat.train.mean)

## check1=mc.pwbart(x.train, post[[6]]$treedraws, post[[6]]$mu, mc.cores=C)
## plot(apply(check1, 2, mean), post[[6]]$yhat.train.mean)

## check2=mc.pwbart(x.train, post[[6]]$treedraws, post[[6]]$mu, mc.cores=C,
##                  dodraws=FALSE)
## plot(check2, post[[6]]$yhat.train.mean)
