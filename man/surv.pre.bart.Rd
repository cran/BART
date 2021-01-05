\name{surv.pre.bart}
\title{Data construction for survival analysis with BART}
\alias{surv.pre.bart}
\description{
Survival data contained in \eqn{(t,\delta, x)} must be translated to data
suitable for the BART survival analysis model; see \code{surv.bart} for more details.
}
\usage{
surv.pre.bart( times, delta, x.train=NULL, x.test=NULL,
               K=NULL, events=NULL, ztimes=NULL, zdelta=NULL )
}
\arguments{
   \item{times}{
     The time of event or right-censoring.\cr
   }

   \item{delta}{
     The event indicator: 1 is an event while 0 is censored.\cr
   }

   \item{x.train}{
   Explanatory variables for training (in sample) data.\cr
   If provided, must be a matrix
   with (as usual) rows corresponding to observations and columns to variables.\cr
   }

   \item{x.test}{
   Explanatory variables for test (out of sample) data.\cr
   If provided, must be a matrix and have the same structure as x.train.\cr
   }

   \item{K}{
     If provided, then coarsen \code{times} per the quantiles
     \eqn{1/K, 2/K, ..., K/K}.
   }
  
   \item{events}{
     If provided, then use for the grid of time points.
   }

   \item{ztimes}{
     If provided, then these columns of \code{x.train} (and
     \code{x.test} if any) are the times for time-dependent covariates.
     They will be transformed into time-dependent covariate sojourn times.
   }
   
   \item{zdelta}{
     If provided, then these columns of \code{x.train} (and
     \code{x.test} if any) are the delta for time-dependent covariates.
     They will be transformed into time-dependent covariate binary events.
   }
}

\value{
   \code{surv.pre.bart} returns a list.
Besides the items listed below, the list has
 a \code{times} component giving the unique times and \code{K} which is the number of
 unique times.
   \item{y.train}{A vector of binary responses.}
   \item{tx.train}{A matrix with rows consisting of time and the covariates of the training data.}
   \item{tx.test}{A matrix with rows consisting of time and the
   covariates of the test data, if any.}
}
%% \references{
%%   Sparapani, R., Logan, B., McCulloch, R., and Laud, P. (2016)
%%   Nonparametric survival analysis using Bayesian Additive Regression Trees (BART).
%%    \emph{Statistics in Medicine}, \bold{16:2741-53} <doi:10.1002/sim.6893>.
%% }
%% \author{
%% Rodney Sparapani: \email{rsparapa@mcw.edu}
%% }
\seealso{
\code{\link{surv.bart}}
}
\examples{

## load the advanced lung cancer example
data(lung)

group <- -which(is.na(lung[ , 7])) ## remove missing row for ph.karno
times <- lung[group, 2]   ##lung$time
delta <- lung[group, 3]-1 ##lung$status: 1=censored, 2=dead
                          ##delta: 0=censored, 1=dead

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

pre <- surv.pre.bart(times=times, delta=delta, x.train=x.train, x.test=x.test)
str(pre)
    
}
\keyword{data construction}
