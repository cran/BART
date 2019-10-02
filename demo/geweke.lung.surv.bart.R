
library(survival)
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
times <- ceiling(times/30) ## months
##times <- ceiling(times/7)  ## weeks

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
post <- mc.surv.bart(x.train=x.train, times=times, delta=delta, x.test=x.train,
                     mc.cores=8, seed=99)

K <- post$K
## select 10 lung cancer patients uniformly spread out over the data set
h <- seq(1, N*K, floor(N/10)*K)

for(i in h) {
    post.mcmc <- post$yhat.test[ , (i-1)+1:K]
    z <- gewekediag(post.mcmc)$z
    y <- max(c(4, abs(z)))

    ## plot the z scores vs. time for each patient
    if(i==1) plot(post$times, z, ylim=c(-y, y), type='l',
                  xlab='t', ylab='z')
    else lines(post$times, z, type='l')
}
## add two-sided alpha=0.05 critical value lines
lines(post$times, rep(-1.96, K), type='l', lty=2)
lines(post$times, rep( 1.96, K), type='l', lty=2)

