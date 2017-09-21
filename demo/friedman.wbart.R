
library(BART)

##simulate from Friedman's five-dimensional test function
##Friedman JH. Multivariate adaptive regression splines
##(with discussion and a rejoinder by the author).
##Annals of Statistics 1991; 19:1-67.

f = function(x) #only the first 5 matter
10*sin(pi*x[ , 1]*x[ , 2]) + 20*(x[ , 3]-.5)^2+10*x[ , 4]+5*x[ , 5]

sigma = 1.0  #y = f(x) + sigma*z where z~N(0, 1)
n = 10000    #number of observations
k = 500      #number of covariates

set.seed(12)
x.train=matrix(runif(n*k), n, k)
Ey.train = f(x.train)
y.train=Ey.train+sigma*rnorm(n)

set.seed(21)
m = 10000
x.test=matrix(runif(m*k), m, k)
Ey.test = f(x.test)
y.test=Ey.test+sigma*rnorm(m)

ndpost = 1000
nskip = 100
C = 8

##run BART with 1 core
num <- proc.time()
set.seed(99)
train = wbart(x.train, y.train, ndpost=ndpost, nskip=nskip)
num <- proc.time()-num
print(num[3])

##run BART with C cores in parallel
den <- proc.time()
mc.train = mc.wbart(x.train, y.train, mc.cores=C,
                  seed=99, ndpost=ndpost, nskip=nskip)
den <- proc.time()-den
print(den[3])

print(c("realized gain:"=num[3]/den[3]))

##Amdahl's Law: theoretical maximum gain
b <- nskip/(ndpost+nskip)
print(c("Amdahl's Law:"=1/((1-b)/C+b)))

par(mfrow=c(2, 2))
y.min <- min(y.train)
y.max <- max(y.train)
plot(y.train, train$yhat.train.mean, asp=1, pch='.',
     xlim=c(y.min, y.max), ylab='wbart train', xlab='true train')
lines(c(y.min, y.max), c(y.min, y.max), type='l')
plot(y.train, mc.train$yhat.train.mean, asp=1, pch='.',
     xlim=c(y.min, y.max), ylab='mc.wbart train', xlab='true train')
lines(c(y.min, y.max), c(y.min, y.max), type='l')

##run predict with 1 core
num <- proc.time()
test = predict(train, x.test)
num <- proc.time()-num
print(num[3])

##run predict with C cores in parallel
den <- proc.time()
mc.test = predict(train, x.test, mc.cores=C)
den <- proc.time()-den
print(den[3])

print(c("realized gain:"=num[3]/den[3]))

y.min <- min(y.test)
y.max <- max(y.test)
yhat.test.mean <- apply(test, 2, mean)
plot(y.test, yhat.test.mean, asp=1, pch='.',
     xlim=c(y.min, y.max), ylab='wbart test', xlab='true test')
lines(c(y.min, y.max), c(y.min, y.max), type='l')
mc.yhat.test.mean <- apply(mc.test, 2, mean)
plot(y.test, mc.yhat.test.mean, asp=1, pch='.',
     xlim=c(y.min, y.max), ylab='mc.wbart test', xlab='true test')
lines(c(y.min, y.max), c(y.min, y.max), type='l')
par(mfrow=c(1, 1))
