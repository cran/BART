\name{mc.surv.pwbart}
\title{Predicting new observations with a previously fitted BART model}
\alias{surv.pwbart}
\alias{mc.surv.pwbart}
\alias{recur.pwbart}
\alias{mc.recur.pwbart}
\description{
BART is a Bayesian \dQuote{sum-of-trees} model.\cr
For a numeric response \eqn{y}, we have
\eqn{y = f(x) + \epsilon}{y = f(x) + e},
where \eqn{\epsilon \sim N(0,\sigma^2)}{e ~ N(0,sigma^2)}.\cr

\eqn{f} is the sum of many tree models.
The goal is to have very flexible inference for the uknown
function \eqn{f}.

In the spirit of \dQuote{ensemble models},
each tree is constrained by a prior to be a weak learner
so that it contributes a
small amount to the overall fit.
}
\usage{
surv.pwbart(
                x.test,
                treedraws,
                binaryOffset=0,
                mc.cores=1L,
                type='pbart',
                transposed=FALSE, nice=19L
              )

mc.surv.pwbart(
                x.test,
                treedraws,
                binaryOffset=0,
                mc.cores=2L,
                type='pbart',
                transposed=FALSE, nice=19L
              )

mc.recur.pwbart(
                x.test,
                treedraws,
                binaryOffset=0,
                mc.cores=2L, 
                type='pbart',
                transposed=FALSE, nice=19L
               )

}
\arguments{
   \item{x.test}{
   Matrix of covariates to predict \eqn{y} for.
   }

   \item{binaryOffset}{
   Mean to add on to \eqn{y} prediction.
   }

   \item{treedraws}{
     \code{$treedraws} returned from \code{surv.bart},
     \code{mc.surv.bart}, \code{recur.bart} or \code{mc.recur.bart}.
   }

   \item{mc.cores}{
   Number of threads to utilize.
 }
 
   \item{type}{ Whether to employ Albert-Chib, \code{'pbart'}, or
     Holmes-Held, \code{'lbart'}. }

   \item{transposed}{
   When running \code{pwbart} or \code{mc.pwbart} in parallel, it is more memory-efficient
   to transpose \code{x.test} prior to calling the internal versions of these functions.
   }

   \item{nice}{
     Set the job niceness.  The default
     niceness is 19: niceness goes from 0 (highest) to 19 (lowest).
   }
}
\details{
   BART is an Bayesian MCMC method.
   At each MCMC interation, we produce a draw from the joint posterior
   \eqn{(f,\sigma) | (x,y)}{(f,sigma) \| (x,y)} in the numeric \eqn{y} case
   and just \eqn{f} in the binary \eqn{y} case.

   Thus, unlike a lot of other modelling methods in R, we do not produce a single model object
   from which fits and summaries may be extracted.  The output consists of values
   \eqn{f^*(x)}{f*(x)} (and \eqn{\sigma^*}{sigma*} in the numeric case) where * denotes a particular draw.
   The \eqn{x} is either a row from the training data (x.train) or the test data (x.test).
}
\value{
Returns an object of type \code{survbart} which is essentially a list with components:

   \item{yhat.test}{
   A matrix with ndpost rows and nrow(x.test) columns.
   Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{f}
   and each column corresponds to a row of x.train.
   The \eqn{(i,j)} value is \eqn{f^*(x)}{f*(x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
   and the \eqn{j^{th}}{j\^th} row of x.train.\cr
   Burn-in is dropped.
   }

   %% \item{surv.test}{test data fits for survival probability.}
   %% \item{surv.test.mean}{mean of \code{surv.test} over the posterior samples.}

   \item{surv.test}{test data fits for survival probability: not
   available for \code{mc.recur.pwbart}.}
   \item{surv.test.mean}{mean of \code{surv.test} over the posterior samples: not
   available for \code{mc.recur.pwbart}.}

   \item{haz.test}{test data fits for hazard: available for
   \code{mc.recur.pwbart} only.}
   \item{haz.test.mean}{mean of \code{haz.test} over the posterior samples:
   available for \code{mc.recur.pwbart} only.}

   \item{cum.test}{test data fits for cumulative hazard: available for
   \code{mc.recur.pwbart} only.}
   \item{cum.test.mean}{mean of \code{cum.test} over the posterior samples:
   available for \code{mc.recur.pwbart} only.}
}
%% \references{
%%   Sparapani, R., Logan, B., McCulloch, R., and Laud, P. (2016)
%%   Nonparametric survival analysis using Bayesian Additive Regression Trees (BART).
%%    \emph{Statistics in Medicine}, \bold{16:2741-53} <doi:10.1002/sim.6893>.
%% }
%% \author{
%% Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
%% Rodney Sparapani: \email{rsparapa@mcw.edu}.
%% }
\seealso{
\code{\link{pwbart}}
}
\examples{

## load the advanced lung cancer example
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

## this x.test is relatively small, but often you will want to
## predict for a large x.test matrix which may cause problems
## due to consumption of RAM so we can predict separately

## mcparallel/mccollect do not exist on windows
if(.Platform$OS.type=='unix') {
##test BART with token run to ensure installation works
    set.seed(99)
    post <- surv.bart(x.train=x.train, times=times, delta=delta, nskip=5, ndpost=5, keepevery=1)

    pre <- surv.pre.bart(x.train=x.train, times=times, delta=delta, x.test=x.test)

    pred <- mc.surv.pwbart(pre$tx.test, post$treedraws, post$binaryOffset)
}

\dontrun{
## run one long MCMC chain in one process
set.seed(99)
post <- surv.bart(x.train=x.train, times=times, delta=delta)

## run "mc.cores" number of shorter MCMC chains in parallel processes
## post <- mc.surv.bart(x.train=x.train, times=times, delta=delta,
##                      mc.cores=8, seed=99)

pre <- surv.pre.bart(x.train=x.train, times=times, delta=delta, x.test=x.test)

pred <- surv.pwbart(pre$tx.test, post$treedraws, post$binaryOffset)

## let's look at some survival curves
## first, a younger group with a healthier KPS
## age 50 with KPS=90: males and females
## males: row 17, females: row 23
x.test[c(17, 23), ]

low.risk.males <- 16*post$K+1:post$K ## K=unique times including censoring
low.risk.females <- 22*post$K+1:post$K

plot(post$times, pred$surv.test.mean[low.risk.males], type='s', col='blue',
     main='Age 50 with KPS=90', xlab='t', ylab='S(t)', ylim=c(0, 1))
points(post$times, pred$surv.test.mean[low.risk.females], type='s', col='red')

}
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
