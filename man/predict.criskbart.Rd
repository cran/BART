\name{predict.criskbart}
\title{Predicting new observations with a previously fitted BART model}
\alias{predict.criskbart}
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
\method{predict}{criskbart}(object, newdata, newdata2, mc.cores=1, openmp=(mc.cores.openmp()>0), ...)
}
\arguments{
   \item{object}{
     \code{object} returned from previous BART fit with \code{crisk.bart}
     or \code{mc.crisk.bart}.
   }

   \item{newdata}{
   Matrix of covariates to predict the distribution of \eqn{t1}.
   }

   \item{newdata2}{
   Matrix of covariates to predict the distribution of \eqn{t2}.
   }

   \item{mc.cores}{
   Number of threads to utilize.
   }

   \item{openmp}{
     Logical value dictating whether OpenMP is utilized for parallel
     processing.  Of course, this depends on whether OpenMP is available
     on your system which, by default, is verified with \code{mc.cores.openmp}.
   }

   \item{...}{
     Other arguments which will be passed on to \code{pwbart}.
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

  Returns an object of type \code{criskbart} with predictions
  corresponding to \code{newdata} and \code{newdata2}.
}
\references{
Chipman, H., George, E., and McCulloch R. (2010)
   Bayesian Additive Regression Trees.
   \emph{The Annals of Applied Statistics}, \bold{4,1}, 266-298 <doi:10.1214/09-AOAS285>.

Chipman, H., George, E., and McCulloch R. (2006)
   Bayesian Ensemble Learning.
   Advances in Neural Information Processing Systems 19,
   Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.

Friedman, J.H. (1991)
   Multivariate adaptive regression splines.
        \emph{The Annals of Statistics}, \bold{19}, 1--67.
}
\author{
Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
Rodney Sparapani: \email{rsparapa@mcw.edu}.
}
\seealso{
\code{\link{crisk.bart}}, \code{\link{mc.crisk.bart}}, \code{\link{mc.crisk.pwbart}}, \code{\link{mc.cores.openmp}}
}
\examples{
##simulate data (example from Friedman MARS paper)
f = function(x){
10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma = 1.0  #y = f(x) + sigma*z , z~N(0,1)
n = 100      #number of observations
set.seed(99)
x=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
y=f(x)

##test BART with token run to ensure installation works
set.seed(99)
post = wbart(x,y,nskip=5,ndpost=5)
x.test = matrix(runif(500*10),500,10)

\dontrun{
##run BART
set.seed(99)
post = wbart(x,y)
x.test = matrix(runif(500*10),500,10)
pred = predict(post, x.test, mu=mean(y))

plot(apply(pred, 2, mean), f(x.test))

}
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}