
library(BART)

## load survival package for the advanced lung cancer example
data(lung)

N <- length(lung$status)

table(lung$ph.karno, lung$pat.karno)

h <- which(is.na(lung$ph.karno))
lung$ph.karno[h] <- lung$pat.karno[h]

times <- lung$time
delta <- lung$status-1 ##lung$status: 1=censored, 2=dead
##delta: 0=censored, 1=dead

## this study reports time in days rather than weeks or months 
## coarsening from days to weeks or months will reduce the computational burden
##times <- ceiling(times/30)
times <- ceiling(times/7)  ## weeks

table(times)
table(delta)

## matrix of observed covariates
x.train <- cbind(lung$age, lung$sex, lung$ph.karno)

## lung$age:        Age in years
## lung$sex:        Male=1 Female=2
## lung$ph.karno:   Karnofsky performance score (dead=0:normal=100:by=10)
##                  rated by physician

dimnames(x.train)[[2]] <- c('age(yr)', 'M(1):F(2)', 'ph.karno(0:100:10)')

## run one long MCMC chain in one process
## set.seed(99)
## post <- surv.bart(x.train=x.train, times=times, delta=delta, x.test=x.test)

## in the interest of time, consider speeding it up by parallel processing
## run "mc.cores" number of shorter MCMC chains in parallel processes
post <- mc.surv.bart(x.train=x.train, times=times, delta=delta, 
                     mc.cores=8, seed=99)

summary(x.train[ , 1])
table(x.train[ , 2])
table(x.train[ , 3])

pre <- surv.pre.bart(times=times, delta=delta, x.train=x.train,
                     x.test=x.train)

pre$tx.test[ , 3] <- 1 ## N.B. col 3, not 2, since time is in col 1
male <- predict(post, newdata=pre$tx.test, mc.cores=8)

pre$tx.test[ , 3] <- 2
female <- predict(post, newdata=pre$tx.test, mc.cores=8)

K <- pre$K

pd.m <- matrix(nrow=1000, ncol=K)
pd.f <- matrix(nrow=1000, ncol=K)

for(j in 1:K) {
    h <- seq(j, N*K, by=K)

    pd.m[ , j] <- apply(male$surv.test[ , h], 1, mean)
    pd.f[ , j] <- apply(female$surv.test[ , h], 1, mean)
}

pd.m.mu <- c(1, apply(pd.m, 2, mean))
pd.m.025 <- c(1, apply(pd.m, 2, quantile, probs=0.025))
pd.m.975 <- c(1, apply(pd.m, 2, quantile, probs=0.975))

pd.f.mu <- c(1, apply(pd.f, 2, mean))
pd.f.025 <- c(1, apply(pd.f, 2, quantile, probs=0.025))
pd.f.975 <- c(1, apply(pd.f, 2, quantile, probs=0.975))

plot(c(0, pre$times), pd.m.mu, col='blue', type='s',
     ylim=0:1, ylab='S(t, x)', xlab='t (weeks)')
lines(c(0, pre$times), pd.m.025, col='blue', type='s', lty=2)
lines(c(0, pre$times), pd.m.975, col='blue', type='s', lty=2)
lines(c(0, pre$times), pd.f.mu, col='red', type='s')
lines(c(0, pre$times), pd.f.025, col='red', type='s', lty=2)
lines(c(0, pre$times), pd.f.975, col='red', type='s', lty=2)

