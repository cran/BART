
library(BART)

B <- getOption('mc.cores', 1)
figures = getOption('figures', default='NONE')

## load survival package for the advanced lung cancer example
data(lung)

N <- length(lung$status)

table(lung$ph.karno, lung$pat.karno)

## if physician's KPS unavailable, then use the patient's
h <- which(is.na(lung$ph.karno))
lung$ph.karno[h] <- lung$pat.karno[h]

times <- lung$time
delta <- lung$status-1 ##lung$status: 1=censored, 2=dead
##delta: 0=censored, 1=dead

## this study reports time in days rather than weeks or months
## coarsening from days to weeks or months will reduce the computational burden
##times <- ceiling(times/30)
times <- ceiling(times/7)  ## weeks

##table(times)
table(delta)

## matrix of observed covariates
x.train <- cbind(lung$sex, lung$age, lung$ph.karno)

## lung$sex:        Male=1 Female=2
## lung$age:        Age in years
## lung$ph.karno:   Karnofsky performance score (dead=0:normal=100:by=10)
##                  rated by physician

dimnames(x.train)[[2]] <- c('M(1):F(2)', 'age(39:82)', 'ph.karno(50:100:10)')

table(x.train[ , 1])
summary(x.train[ , 2])
table(x.train[ , 3])

## run one long MCMC chain in one process
## set.seed(99)
## post <- surv.bart(x.train=x.train, times=times, delta=delta, x.test=x.test)

## in the interest of time, consider speeding it up by parallel processing
## run "mc.cores" number of shorter MCMC chains in parallel processes
post <- mc.surv.bart(x.train=x.train, times=times, delta=delta,
                     mc.cores=B, seed=99, K=100)

pre <- surv.pre.bart(times=times, delta=delta, x.train=x.train,
                     x.test=x.train, K=100)

K <- pre$K
M <- post$ndpost
NK <- N*K

pre$tx.test <- rbind(pre$tx.test, pre$tx.test)
pre$tx.test[ , 2] <- c(rep(1, N*K), rep(2, N*K))
## sex pushed to col 2, since time is always in col 1

pred <- predict(post, newdata=pre$tx.test, mc.cores=B)

for(i in seq(1, N, by=5)) {
##for(i in 1:N) {
    h=(i-1)*K+1:K
    if(i==1)
        plot(c(0, pre$times), c(1, pred$surv.test.mean[h]),
             type='s', col=4, lty=2,
             ylim=0:1, ylab='S(t, x)', xlab='t (weeks)',)
    else lines(c(0, pre$times), c(1, pred$surv.test.mean[h]),
               type='s', col=4, lty=2)
    lines(c(0, pre$times), c(1, pred$surv.test.mean[h+NK]),
          type='s', col=2, lty=3)
}

if(figures!='NONE')
    dev.copy2pdf(file=paste(figures, 'lung-ice.pdf', sep='/'))
