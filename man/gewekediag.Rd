\name{gewekediag}
\alias{gewekediag}
\title{Geweke's convergence diagnostic}

\usage{gewekediag(x, frac1=0.1, frac2=0.5)}

\arguments{
  \item{x}{Matrix of MCMC chains: the rows are the samples and
  the columns are different "parameters".  For BART, generally, the
  columns are estimates of \eqn{f}.  For \code{pbart}, they are
  different subjects.  For \code{surv.bart}, they are different subjects
  at a grid of times.}
  \item{frac1}{fraction to use from beginning of chain}
  \item{frac2}{fraction to use from end of chain}
}

\value{
Z-scores for a test of equality of means
between the first and last parts of the chain. A separate
statistic is calculated for each variable in each chain.
}

\description{
Geweke (1992) proposed a convergence diagnostic for Markov chains
based on a test for equality of the means of the first and last part
of a Markov chain (by default the first 10\% and the last 50\%).  If the
samples are drawn from the stationary distribution of the chain, the two
means are equal and Geweke's statistic has an asymptotically standard
normal distribution.

The test statistic is a standard Z-score: the difference between the
two sample means divided by its estimated standard error.  The standard
error is estimated from the spectral density at zero and so takes into
account any autocorrelation.

The Z-score is calculated under the assumption that the two parts of
the chain are asymptotically independent, which requires that the sum
of \code{frac1} and \code{frac2} be strictly less than 1.

Adapted from the \code{geweke.diag} function of
the coda package which passes \code{mcmc} objects as arguments
rather than matrices.
}

\seealso{
   \code{\link{spectrum0ar}}.
}

\references{
  Geweke J. (1992)
  Evaluating the Accuracy of Sampling-Based Approaches to the Calculation of Posterior Moments.
  In JM Bernado, JO Berger, AP Dawid, AFM Smith (eds.), Bayesian Statistics 4, pp. 169-193.
  Oxford University Press, Oxford.
  
  Plummer M., Best N., Cowles K. and Vines K. (2006)
  CODA: Convergence Diagnosis and Output Analysis for MCMC.
  R News, vol 6, 7-11.
}
\examples{

## load survival package for the advanced lung cancer example
data(lung)

group <- -which(is.na(lung[ , 7])) ## remove missing row for ph.karno
times <- lung[group, 2]   ##lung$time
delta <- lung[group, 3]-1 ##lung$status: 1=censored, 2=dead
                          ##delta: 0=censored, 1=dead

## this study reports time in days rather than months like other studies
## coarsening from days to months will reduce the computational burden
times <- ceiling(times/30)

summary(times)
table(delta)

x.train <- as.matrix(lung[group, c(4, 5, 7)]) ## matrix of observed covariates

## lung$age:        Age in years
## lung$sex:        Male=1 Female=2
## lung$ph.karno:   Karnofsky performance score (dead=0:normal=100:by=10)
##                  rated by physician

dimnames(x.train)[[2]] <- c('age(yr)', 'M(1):F(2)', 'ph.karno(0:100:10)')

summary(x.train[ , 1])
table(x.train[ , 2])
table(x.train[ , 3])

x.test <- matrix(nrow=84, ncol=3) ## matrix of covariate scenarios

dimnames(x.test)[[2]] <- dimnames(x.train)[[2]]

i <- 1

for(age in 5*(9:15)) for(sex in 1:2) for(ph.karno in 10*(5:10)) {
    x.test[i, ] <- c(age, sex, ph.karno)
    i <- i+1
}

\dontrun{
    set.seed(99)
    post <- surv.bart(x.train=x.train, times=times, delta=delta, x.test=x.test)
    ## in the interest of time, consider speeding it up by parallel processing
    ## run "mc.cores" number of shorter MCMC chains in parallel processes
    ## post <- mc.surv.bart(x.train=x.train, times=times, delta=delta,
    ##                      x.test=x.test, mc.cores=8, seed=99)

    N <- nrow(x.test)

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

}

}

\keyword{convergence diagnostics}
