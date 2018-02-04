
library(BART)

##simulate from Friedman's five-dimensional test function
##Friedman JH. Multivariate adaptive regression splines
##(with discussion and a rejoinder by the author).
##Annals of Statistics 1991; 19:1-67.

f = function(x) #only the first 5 matter
    sin(pi*x[ , 1]*x[ , 2]) + 2*(x[ , 3]-.5)^2+x[ , 4]+0.5*x[ , 5]-1.5

sigma = 1.0  #y = f(x) + sigma*z where z~N(0, 1)
P = 100      #number of covariates
C = 8
thin <- c(10, 50, 250)

par(mfrow=c(3, 1))

post <- as.list(1:3)

for(i in 1:3) {
    N <- 10^(1+i)
    set.seed(12)
    x.train=matrix(runif(N*P), N, P)
    Ey.train = f(x.train)
    y.train=((Ey.train+sigma*rnorm(N))>0)*1

    post[[i]] = mc.lbart(x.train, y.train, mc.cores=C,
                         keepevery=1, sparse=TRUE, seed=99)

    plot(post[[i]]$varprob.mean, col=c(rep(2, 5), rep(1, P-5)),
         main=paste0('N:', N, ', P:', P, ', thin:', thin[i]),
         ylab='Selection Probability', ylim=c(0, 0.2),
         pch=1+45*(post[[i]]$varprob.mean <= 1/P))
    lines(c(0, 100), c(1/P, 1/P))

    table(1+45*(post[[i]]$varprob.mean <= 1/P))
}

par(mfrow=c(1, 1))

##dev.copy2pdf(file='sparse-lbart.pdf')
