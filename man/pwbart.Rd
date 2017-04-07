\name{pwbart}
\title{Predicting new observations with a previously fitted BART model}
\alias{pwbart}
\alias{mc.pwbart}
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
pwbart( x.test, treedraws, mu=0, mc.cores=1L, transposed=FALSE,
        nice=19L ## mc.pwbart only
      )

mc.pwbart( x.test, treedraws, mu=0, mc.cores=2L, transposed=FALSE,
           nice=19L ## mc.pwbart only
         )
}
\arguments{
   \item{x.test}{
   Matrix of covariates to predict \eqn{y} for.
   }

   \item{treedraws}{
     \code{$treedraws} returned from \code{wbart} or \code{pbart}.
   }

   \item{mu}{
   Mean to add on to \eqn{y} prediction.
   }

   \item{mc.cores}{
   Number of threads to utilize.
   }

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
   The \code{plot} method sets mfrow to c(1,2) and makes two plots.\cr
   The first plot is the sequence of kept draws of \eqn{\sigma}{sigma}
   including the burn-in draws.  Initially these draws will decline as BART finds fit
   and then level off when the MCMC has burnt in.\cr
   The second plot has \eqn{y} on the horizontal axis and posterior intervals for
   the corresponding \eqn{f(x)} on the vertical axis.


   \code{bart} returns a list assigned class \sQuote{bart}.
   In the numeric \eqn{y} case, the list has components:

   \item{yhat.train}{
   A matrix with (ndpost/keepevery) rows and nrow(x.train) columns.
   Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{f}
   and each column corresponds to a row of x.train.
   The \eqn{(i,j)} value is \eqn{f^*(x)}{f*(x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
   and the \eqn{j^{th}}{j\^th} row of x.train.\cr
   Burn-in is dropped.
   }

   \item{yhat.test}{Same as yhat.train but now the x's are the rows of the test data.}
   \item{yhat.train.mean}{train data fits = mean of yhat.train columns.}
   \item{yhat.test.mean}{test data fits = mean of yhat.test columns.}
   \item{sigma}{post burn in draws of sigma, length = ndpost/keepevery.}
   \item{first.sigma}{burn-in draws of sigma.}
   \item{varcount}{a matrix with (ndpost/keepevery) rows and nrow(x.train) columns.
   Each row is for a draw. For each variable (corresponding to the columns),
   the total count of the number of times
   that variable is used in a tree decision rule (over all trees) is given.}

   \item{sigest}{
   The rough error standard deviation (\eqn{\sigma}{sigma}) used in the prior.
   }

   \item{y}{
   The input dependent vector of values for the dependent variable.\cr
   This is used in plot.bart.
   }

   In the binary \eqn{y} case, the returned list has the components
   yhat.train, yhat.test, and varcount as above.  In addition the list
   has a binaryOffset component giving the value used.

   Note that in the binary \eqn{y}, case yhat.train and yhat.test are
   \eqn{f(x)} + binaryOffset.  If you want draws of the probability
   \eqn{P(Y=1 | x)} you need to apply the normal cdf (\code{pnorm})
   to these values.
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
\code{\link{wbart}}
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
pred = pwbart(post$treedraws, x.test, mu=mean(y))

plot(apply(pred, 2, mean), f(x.test))

}
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
