\name{mc.lbart}
\title{Logit BART for dichotomous outcomes with Logistic latents and parallel computation}
\alias{mc.lbart}
\description{
BART is a Bayesian \dQuote{sum-of-trees} model.\cr
For numeric response \eqn{y}, we have
\eqn{y = f(x) + \epsilon}{y = f(x) + e},
where \eqn{\epsilon \sim Log(0, 1)}{e ~ Log(0, 1)}.\cr
For a binary response \eqn{y}, \eqn{P(Y=1 | x) = F(f(x))}, where \eqn{F}
denotes the standard Logistic CDF (logit link).

In both cases, \eqn{f} is the sum of many tree models.
The goal is to have very flexible inference for the uknown
function \eqn{f}.

In the spirit of \dQuote{ensemble models},
each tree is constrained by a prior to be a weak learner
so that it contributes a
small amount to the overall fit.
}
\usage{
mc.lbart(
   x.train, y.train, x.test=matrix(0.0,0,0),
   sparse=FALSE, a=0.5, b=1, augment=FALSE, rho=NULL,
   xinfo=matrix(0.0,0,0), usequants=FALSE,
   cont=FALSE, rm.const=TRUE, tau.interval=0.95,
   k=2.0, power=2.0, base=.95,
   binaryOffset=NULL,
   ntree=50L, numcut=100L,
   ndpost=1000L, nskip=100L,
   keepevery=1L, printevery=100,
   keeptrainfits=TRUE, transposed=FALSE,
   %treesaslists = FALSE,
   mc.cores = 2L, nice = 19L,
   seed = 99L
)
}
\arguments{
   \item{x.train}{
   Explanatory variables for training (in sample) data.\cr
   May be a matrix or a data frame,
   with (as usual) rows corresponding to observations and columns to variables.\cr
   If a variable is a factor in a data frame, it is replaced with dummies.
   Note that q dummies are created if q>2 and
   one dummy is created if q=2, where q is the number of levels of the factor.
   \code{lbart} will generate draws of \eqn{f(x)} for each \eqn{x}
   which is a row of x.train.
   }

   \item{y.train}{
   Dependent variable for training (in sample) data.\cr
   If y is numeric a continous response model is fit (normal errors).\cr
   If y is a factor (or just has values 0 and 1) then a binary response model
   with a logit link is fit.
   }

   \item{x.test}{
   Explanatory variables for test (out of sample) data.\cr
   Should have same structure as x.train.\cr
   \code{lbart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of x.test.
   }

   \item{sparse}{Whether to perform variable selection based on a
     sparse Dirichlet prior rather than simply uniform; see Linero 2016.}
   \item{a}{Sparse parameter for \eqn{Beta(a, b)} prior:
     \eqn{0.5<=a<=1} where lower values inducing more sparsity.}
   \item{b}{Sparse parameter for \eqn{Beta(a, b)} prior; typically,
     \eqn{b=1}.}
   \item{rho}{Sparse parameter: typically \eqn{rho=p} where \eqn{p} is the
     number of covariates under consideration.}
   \item{augment}{Whether data augmentation is to be performed in sparse
     variable selection.}
   
   \item{xinfo}{ You can provide the cutpoints to BART or let BART
     choose them for you.  To provide them, use the \code{xinfo}
     argument to specify a list (matrix) where the items (rows) are the
     covariates and the contents of the items (columns) are the
     cutpoints.  }

   \item{usequants}{ If \code{usequants=FALSE}, then the
    cutpoints in \code{xinfo} are generated uniformly; otherwise,
    if \code{TRUE}, uniform quantiles are used for the cutpoints. }
   
   \item{cont}{ Whether or not to assume all variables are continuous.}

   \item{rm.const}{ Whether or not to remove constant variables.}  

   \item{tau.interval}{
     The width of the interval to scale the variance for the terminal
     leaf values.}
   
   \item{k}{
   For numeric y,
   k is the number of prior standard deviations \eqn{E(Y|x) = f(x)} is away from +/-.5.
   The response (y.train) is internally scaled to range from -.5 to .5.
   For binary y,
   k is the number of prior standard deviations \eqn{f(x)} is away from +/-3.
   In both cases, the bigger k is, the more conservative the fitting will be.
   }

   \item{power}{
   Power parameter for tree prior.
   }

   \item{base}{
   Base parameter for tree prior.
   }

   \item{binaryOffset}{
   Used for binary \eqn{y}.\cr
   The model is \eqn{P(Y=1 | x) = F(f(x) + binaryOffset)}.\cr
   %% The idea is that \eqn{f} is shrunk towards 0, so the offset allows you to shrink towards
   %% a probability other than .5.
   }

   \item{ntree}{
   The number of trees in the sum.
   }

   \item{numcut}{
   The number of possible values of c (see usequants).
   If a single number if given, this is used for all variables.
   Otherwise a vector with length equal to ncol(x.train) is required,
   where the \eqn{i^{th}}{i^th} element gives the number of c used for
   the \eqn{i^{th}}{i^th} variable in x.train.
   If usequants is false, numcut equally spaced cutoffs
   are used covering the range of values in the corresponding
   column of x.train.  If usequants is true, then  min(numcut, the number of unique values in the
   corresponding columns of x.train - 1) c values are used.
   }

   \item{ndpost}{
   The number of posterior draws returned.
   }

   \item{nskip}{
   Number of MCMC iterations to be treated as burn in.
   }

   \item{keepevery}{
   Every keepevery draw is kept to be returned to the user.
   }

   \item{printevery}{
   As the MCMC runs, a message is printed every printevery draws.
   }

   \item{keeptrainfits}{
   Whether to keep \code{yhat.train} or not.
   }

   \item{transposed}{
   When running \code{lbart} in parallel, it is more memory-efficient
   to transpose \code{x.train} and \code{x.test}, if any, prior to
   calling \code{mc.lbart}.
 }

   %% \item{treesaslists}{
   %%   The default behavior is to return the trees as a text string.  If you
   %%   specify \code{treesaslists=TRUE}, then the trees will be returned
   %%   as a list as well.
   %% }

   \item{seed}{
     Setting the seed required for reproducible MCMC.
   }

   \item{mc.cores}{
     Number of cores to employ in parallel.
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
   \code{mc.lbart} returns an object of type \code{lbart} which is
   essentially a list. % assigned class \sQuote{bart}.

   \item{yhat.train}{
   A matrix with ndpost rows and nrow(x.train) columns.
   Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{f}
   and each column corresponds to a row of x.train.
   The \eqn{(i,j)} value is \eqn{f^*(x)}{f*(x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
   and the \eqn{j^{th}}{j\^th} row of x.train.\cr
   Burn-in is dropped.
   }

   \item{yhat.test}{Same as yhat.train but now the x's are the rows of the test data.}
   \item{yhat.train.mean}{train data fits = mean of yhat.train columns.}
   \item{yhat.test.mean}{test data fits = mean of yhat.test columns.}
   \item{varcount}{a matrix with ndpost rows and nrow(x.train) columns.
   Each row is for a draw. For each variable (corresponding to the columns),
   the total count of the number of times
   that variable is used in a tree decision rule (over all trees) is given.}

   In addition, the list
   has a \code{binaryOffset} giving the value used.

   Note that in the binary \eqn{y}, case yhat.train and yhat.test are
   \eqn{f(x) + binaryOffset}.  If you want draws of the probability
   \eqn{P(Y=1 | x)} you need to apply the Logistic cdf (\code{plogis})
   to these values.
}
%% \references{
%% Chipman, H., George, E., and McCulloch R. (2010)
%%    Bayesian Additive Regression Trees.
%%    \emph{The Annals of Applied Statistics}, \bold{4,1}, 266-298 <doi: 10.1214/09-AOAS285>.

%% Friedman, J.H. (1991)
%%    Multivariate adaptive regression splines.
%%    \emph{The Annals of Statistics}, \bold{19}, 1--67.

%% Gramacy, RB and Polson, NG (2012)
%%    Simulation-based regularized logistic regression.
%%    \emph{Bayesian Analysis}, \bold{7}, 567--590.

%% Holmes, C and Held, L (2006)
%%    Bayesian auxiliary variable models for binary and multinomial regression.
%%    \emph{Bayesian Analysis}, \bold{1}, 145--68.
	
%% Linero, A.R. (2018)
%%   Bayesian regression trees for high dimensional prediction and variable
%%   selection. \emph{JASA}, \bold{113}, 626--36.

%% }
%% \author{
%% Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
%% Rodney Sparapani: \email{rsparapa@mcw.edu}.
%% }
\seealso{
\code{\link{lbart}}
}
\examples{
set.seed(99)
n=5000
x = sort(-2+4*runif(n))
X=matrix(x,ncol=1)
f = function(x) {return((1/2)*x^3)}
FL = function(x) {return(exp(x)/(1+exp(x)))}
pv = FL(f(x))
y = rbinom(n,1,pv)
np=100
xp=-2+4*(1:np)/np
Xp=matrix(xp,ncol=1)

## parallel::mcparallel/mccollect do not exist on windows
## if(.Platform$OS.type=='unix') {
## ##test BART with token run to ensure installation works
##     mf = mc.lbart(X, y, nskip=5, ndpost=5, mc.cores=1, seed=99)
## }

\dontrun{
set.seed(99)
pf = lbart(X,y,Xp)

plot(f(Xp), pf$yhat.test.mean, xlim=c(-4, 4), ylim=c(-4, 4),
     xlab='True f(x)', ylab='BART f(x)')
lines(c(-4, 4), c(-4, 4))

mf = mc.lbart(X,y,Xp, mc.cores=4, seed=99)

plot(f(Xp), mf$yhat.test.mean, xlim=c(-4, 4), ylim=c(-4, 4),
     xlab='True f(x)', ylab='BART f(x)')
lines(c(-4, 4), c(-4, 4))

par(mfrow=c(2,2))

plot(range(xp),range(pf$yhat.test),xlab='x',ylab='f(x)',type='n')
lines(x,f(x),col='blue',lwd=2)
lines(xp,apply(pf$yhat.test,2,mean),col='red')
qpl = apply(pf$yhat.test,2,quantile,probs=c(.025,.975))
lines(xp,qpl[1,],col='green',lty=1)
lines(xp,qpl[2,],col='green',lty=1)
title(main='BART::lbart f(x) with 0.95 intervals')

plot(range(xp),range(mf$yhat.test),xlab='x',ylab='f(x)',type='n')
lines(x,f(x),col='blue',lwd=2)
lines(xp,apply(mf$yhat.test,2,mean),col='red')
qpl = apply(mf$yhat.test,2,quantile,probs=c(.025,.975))
lines(xp,qpl[1,],col='green',lty=1)
lines(xp,qpl[2,],col='green',lty=1)
title(main='BART::mc.lbart f(x) with 0.95 intervals')

plot(pf$yhat.test.mean,apply(mf$yhat.test,2,mean),xlab='BART::lbart',ylab='BART::mc.lbart')
abline(0,1,col='red')
title(main="BART::lbart f(x) vs. BART::mc.lbart f(x)")
}
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
