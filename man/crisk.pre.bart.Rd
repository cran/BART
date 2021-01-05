\name{crisk.pre.bart}
\title{Data construction for competing risks with BART}
\alias{crisk.pre.bart}
\description{
Competing risks contained in \eqn{(t, \delta, x)} must be translated to data
suitable for the BART competing risks model; see \code{crisk.bart} for more details.
}
\usage{
crisk.pre.bart( times, delta, x.train=NULL, x.test=NULL,
                x.train2=x.train, x.test2=x.test, K=NULL )
}
\arguments{
   \item{times}{
     The time of event or right-censoring.\cr
   }

   \item{delta}{
     The event indicator: 1 is a cause 1 event, 2 a cause 2 while 0 is censored.\cr
   }

   \item{x.train}{
   Explanatory variables for training (in sample) data of cause 1.\cr
   If provided, must be a matrix
   with (as usual) rows corresponding to observations and columns to variables.\cr
   }

   \item{x.test}{
   Explanatory variables for test (out of sample) data of cause 1.\cr
   If provided, must be a matrix and have the same structure as x.train.\cr
 }
 
   \item{x.train2}{
   Explanatory variables for training (in sample) data of cause 2.\cr
   If provided, must be a matrix
   with (as usual) rows corresponding to observations and columns to variables.\cr
   }

   \item{x.test2}{
   Explanatory variables for test (out of sample) data of cause 2.\cr
   If provided, must be a matrix and have the same structure as x.train.\cr
 }
 
   \item{K}{
     If provided, then coarsen \code{times} per the quantiles
     \eqn{1/K, 2/K, ..., K/K}.
   }
}

\value{
   \code{surv.pre.bart} returns a list.
Besides the items listed below, the list has
 a \code{times} component giving the unique times and \code{K} which is the number of
 unique times.
   \item{y.train}{A vector of binary responses for cause 1.}
   \item{y.train2}{A vector of binary responses for cause 2.}
   \item{cond}{A vector of indices of \code{y.train} indicating censored subjects.}
   \item{binaryOffset}{The binary offset for \code{y.train}.}
   \item{binaryOffset2}{The binary offset for \code{y.train2}.}
   \item{tx.train}{A matrix with rows consisting of time and the
   covariates of the training data for cause 1.}
   \item{tx.train2}{A matrix with rows consisting of time and the
   covariates of the training data for cause 2.}
   \item{tx.test}{A matrix with rows consisting of time and the
   covariates of the test data, if any, for cause 1.}
   \item{tx.test2}{A matrix with rows consisting of time and the
   covariates of the test data, if any, for cause 2.}
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
\code{\link{crisk.bart}}
}
\examples{

data(transplant)

delta <- (as.numeric(transplant$event)-1)

delta[delta==1] <- 4
delta[delta==2] <- 1
delta[delta>1] <- 2
table(delta, transplant$event)

table(1+floor(transplant$futime/30.5)) ## months
times <- 1+floor(transplant$futime/30.5)

typeO <- 1*(transplant$abo=='O')
typeA <- 1*(transplant$abo=='A')
typeB <- 1*(transplant$abo=='B')
typeAB <- 1*(transplant$abo=='AB')
table(typeA, typeO)

x.train <- cbind(typeO, typeA, typeB, typeAB)

N <- nrow(x.train)

x.test <- x.train

x.test[1:N, 1:4] <- matrix(c(1, 0, 0, 0), nrow=N, ncol=4, byrow=TRUE)

pre <- crisk.pre.bart(x.train=x.train, times=times, delta=delta, x.test=x.test)

}
\keyword{data construction}
