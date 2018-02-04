
library(BART)

##simulate from Friedman's five-dimensional test function
##Friedman JH. Multivariate adaptive regression splines
##(with discussion and a rejoinder by the author).
##Annals of Statistics 1991; 19:1-67.

f = function(x) #only the first 5 matter
    sin(pi*x[ , 1]*x[ , 2]) + 2*(x[ , 3]-.5)^2+x[ , 4]+0.5*x[ , 5]-1.5

sigma = 1.0  # y~Logistic(f(x), sigma)
P = 5        # number of covariates
N = 2500
C = 3        # number of categories
M = 5^P

x <- seq(0, 1, length.out=5)
x.test <- matrix(nrow=M, ncol=P)
dimnames(x.test)[[2]] <- paste0('x', 1:5)
h <- 1
for(x5 in x)
    for(x4 in x)
        for(x3 in x)
            for(x2 in x)
                for(x1 in x) {
                    x.test[h, ] <- c(x1, x2, x3, x4, x5)
                    h <- h+1
                }
Ey.test = f(x.test)

set.seed(12)
x.train=matrix(runif(N*P), N, P)
dimnames(x.train)[[2]] <- paste0('x', 1:5)
Ey.train = f(x.train)
y.train=rlogis(N, Ey.train, sigma)

i <- y.train > -1
j <- y.train > 1

y.train[!i] <- 1
y.train[i] <- 2
y.train[j] <- 3

table(y.train)

## set.seed(99)
## post = mbart(x.train, y.train, x.test)

post = mc.mbart(x.train, y.train, x.test, mc.cores=8, seed=99)

h <- seq(1, C*N, by=C)

print(cor(post$prob.train.mean[h], plogis(-1, Ey.train, sigma))^2)
print(cor(post$prob.train.mean[h+1], plogis(1, Ey.train, sigma)-plogis(-1, Ey.train, sigma))^2)
print(cor(post$prob.train.mean[h+2], plogis(1, Ey.train, sigma, 0))^2)

plot(plogis(-1, Ey.train, sigma), post$prob.train.mean[h], pch='.',
     xlim=0:1, ylim=0:1, xlab='Known P(y=1)', ylab='Est. P(y=1)')
abline(0, 1)

plot(plogis(1, Ey.train, sigma)-plogis(-1, Ey.train, sigma), post$prob.train.mean[h+1], pch='.',
     xlim=0:1, ylim=0:1, xlab='Known P(y=2)', ylab='Est. P(y=2)')
abline(0, 1)

plot(plogis(1, Ey.train, sigma, 0), post$prob.train.mean[h+2], pch='.',
     xlim=0:1, ylim=0:1, xlab='Known P(y=3)', ylab='Est. P(y=3)')
abline(0, 1)

h <- seq(1, C*M, by=C)

print(cor(post$prob.test.mean[h], plogis(-1, Ey.test, sigma))^2)
print(cor(post$prob.test.mean[h+1], plogis(1, Ey.test, sigma)-plogis(-1, Ey.test, sigma))^2)
print(cor(post$prob.test.mean[h+2], plogis(1, Ey.test, sigma, 0))^2)

plot(plogis(-1, Ey.test, sigma), post$prob.test.mean[h], pch='.',
     xlim=0:1, ylim=0:1, xlab='Known P(y=1)', ylab='Est. P(y=1)')
abline(0, 1)

plot(plogis(1, Ey.test, sigma)-plogis(-1, Ey.test, sigma), post$prob.test.mean[h+1], pch='.',
     xlim=0:1, ylim=0:1, xlab='Known P(y=2)', ylab='Est. P(y=2)')
abline(0, 1)

plot(plogis(1, Ey.test, sigma, 0), post$prob.test.mean[h+2], pch='.',
     xlim=0:1, ylim=0:1, xlab='Known P(y=3)', ylab='Est. P(y=3)')
abline(0, 1)
