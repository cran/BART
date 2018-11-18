
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
K = 3        # number of categories

set.seed(12)
x.train=matrix(runif(N*P), N, P)
dimnames(x.train)[[2]] <- paste0('x', 1:5)
Ey.train = f(x.train)
y.train=rnorm(N, Ey.train, sigma)

i <- y.train > -1
j <- y.train > 1

y.train[!i] <- 1
y.train[i] <- 2
y.train[j] <- 3

table(y.train)

## set.seed(99)
## post = mbart(x.train, y.train, x.train)

post = mc.mbart(x.train, y.train, x.train, mc.cores=8, seed=99)

h <- seq(1, K*N, by=K)

print(cor(post$prob.test.mean[h], pnorm(-1, Ey.train, sigma))^2)
print(cor(post$prob.test.mean[h+1], pnorm(1, Ey.train, sigma)-pnorm(-1, Ey.train, sigma))^2)
print(cor(post$prob.test.mean[h+2], pnorm(1, Ey.train, sigma, 0))^2)

par(mfrow=c(2, 2))
plot(pnorm(-1, Ey.train, sigma), post$prob.test.mean[h], pch='.',
     xlim=0:1, ylim=0:1, xlab='Known P(y=1)', ylab='Est. P(y=1)')
abline(0, 1)
plot(pnorm(1, Ey.train, sigma)-pnorm(-1, Ey.train, sigma), post$prob.test.mean[h+1], pch='.',
     xlim=0:1, ylim=0:1, xlab='Known P(y=2)', ylab='Est. P(y=2)')
abline(0, 1)
plot(pnorm(1, Ey.train, sigma, 0), post$prob.test.mean[h+2], pch='.',
     xlim=0:1, ylim=0:1, xlab='Known P(y=3)', ylab='Est. P(y=3)')
abline(0, 1)
par(mfrow=c(1, 1))

