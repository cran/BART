\name{wbart}
\title{BART for continuous outcomes}
\alias{wbart}
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
wbart(
   x.train, y.train, x.test=matrix(0.0,0,0),
   sparse=FALSE, theta=0, omega=1,
   a=0.5, b=1, augment=FALSE, rho=NULL,
   xinfo=matrix(0.0,0,0), usequants=FALSE,
   cont=FALSE, rm.const=TRUE,
   sigest=NA, sigdf=3, sigquant=.90,
   k=2.0, power=2.0, base=.95,
   sigmaf=NA, lambda=NA,
   fmean=mean(y.train), w=rep(1,length(y.train)),
   ntree=200L, numcut=100L,
   ndpost=1000L, nskip=100L, keepevery=1L,
   nkeeptrain=ndpost, nkeeptest=ndpost,
   nkeeptestmean=ndpost, nkeeptreedraws=ndpost,
   printevery=100L, transposed=FALSE %, treesaslists=FALSE
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
   %\code{makeind} is used to generate the dummies.
   \code{wbart} will generate draws of \eqn{f(x)} for each \eqn{x}
   which is a row of \code{x.train}.
   }

   \item{y.train}{
   Continuous dependent variable for training (in sample) data.\cr
   %% If y is numeric a continous response model is fit (normal errors).\cr
   %% If y is a factor (or just has values 0 and 1) then a binary response model
   %% with a probit link is fit.
   }

   \item{x.test}{
   Explanatory variables for test (out of sample) data.\cr
   Should have same structure as x.train.\cr
   \code{wbart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of x.test.
   }

   \item{sparse}{Whether to perform variable selection based on a
     sparse Dirichlet prior rather than simply uniform; see Linero 2016.}
   \item{theta}{Set \eqn{theta} parameter; zero means random.}
   \item{omega}{Set \eqn{omega} parameter; zero means random.}
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
  
   \item{sigest}{
   The prior for the error variance  (\eqn{\sigma^2}{sigma\^2}) is inverted chi-squared
   (the standard conditionally conjugate prior).
   The prior is specified by choosing the degrees of freedom, a rough estimate of the
   corresponding standard deviation and a quantile to put this rough estimate at.
   If sigest=NA then the rough estimate will be the usual least squares estimator.
   Otherwise the supplied value will be used.
   %Not used if y is binary.
   }

   \item{sigdf}{
   Degrees of freedom for error variance prior.
   %Not used if y is binary.
   }

   \item{sigquant}{
   The quantile of the prior that the rough estimate (see sigest) is placed at.
   The closer the quantile is to 1,
   the more aggresive the fit will be as you are putting more prior weight
   on error standard deviations (\eqn{\sigma}{sigma}) less than the rough estimate.
   %Not used if y is binary.
   }

   \item{k}{
   For numeric y,
   k is the number of prior standard deviations \eqn{E(Y|x) = f(x)} is away from +/-.5.
   The response (y.train) is internally scaled to range from -.5 to .5.
   %% For binary y,
   k is the number of prior standard deviations \eqn{f(x)} is away from +/-3.
   %In both cases,
   The bigger k is, the more conservative the fitting will be.
   }

   \item{power}{
   Power parameter for tree prior.
   }

   \item{base}{
   Base parameter for tree prior.
   }

   \item{sigmaf}{
   The SD of f.
   }

   \item{lambda}{
   The scale of the prior for the variance.
   }

   \item{fmean}{
   BART operates on \code{y.train} centered by \code{fmean}.
   }

   \item{w}{
   Vector of weights which multiply the standard deviation.
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

   \item{nkeeptrain}{
   Number of MCMC iterations to be returned for train data.
   }

   \item{nkeeptest}{
   Number of MCMC iterations to be returned for test data.
   }

   \item{nkeeptestmean}{
   Number of MCMC iterations to be returned for test mean.
   }

   \item{nkeeptreedraws}{
   Number of MCMC iterations to be returned for tree draws.
   }

   \item{printevery}{
   As the MCMC runs, a message is printed every printevery draws.
   }

   \item{keepevery}{
   Every keepevery draw is kept to be returned to the user.\cr
   %% A \dQuote{draw} will consist of values of the error standard deviation (\eqn{\sigma}{sigma})
   %% and \eqn{f^*(x)}{f*(x)}
   %% at \eqn{x} = rows from the train(optionally) and test data, where \eqn{f^*}{f*} denotes
   %% the current draw of \eqn{f}.
   }

   \item{transposed}{
   When running \code{wbart} in parallel, it is more memory-efficient
   to transpose \code{x.train} and \code{x.test}, if any, prior to
   calling \code{mc.wbart}.
   }

   %% \item{treesaslists}{
   %%   The default behavior is to return the trees as a text string.  If you
   %%   specify \code{treesaslists=TRUE}, then the trees will be returned
   %%   as a list as well.
   %% }
}
\details{
   BART is an Bayesian MCMC method.
   At each MCMC interation, we produce a draw from the joint posterior
   \eqn{(f,\sigma) | (x,y)}{(f,sigma) \| (x,y)} in the numeric \eqn{y} case.
   %and just \eqn{f} in the binary \eqn{y} case.

   Thus, unlike a lot of other modelling methods in R, we do not produce a single model object
   from which fits and summaries may be extracted.  The output consists of values
   \eqn{f^*(x)}{f*(x)} (and \eqn{\sigma^*}{sigma*} in the numeric case) where * denotes a particular draw.
   The \eqn{x} is either a row from the training data (x.train) or the test data (x.test).
}
\value{
   %% The \code{plot} method sets mfrow to c(1,2) and makes two plots.\cr
   %% The first plot is the sequence of kept draws of \eqn{\sigma}{sigma}
   %% including the burn-in draws.  Initially these draws will decline as BART finds fit
   %% and then level off when the MCMC has burnt in.\cr
   %% The second plot has \eqn{y} on the horizontal axis and posterior intervals for
   %% the corresponding \eqn{f(x)} on the vertical axis.

   \code{wbart} returns an object of type \code{wbart} which is
   essentially a list. % assigned class \sQuote{bart}.
   In the numeric \eqn{y} case, the list has components:

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
   \item{sigma}{post burn in draws of sigma, length = ndpost.}
   \item{first.sigma}{burn-in draws of sigma.}
   \item{varcount}{a matrix with ndpost rows and nrow(x.train) columns.
   Each row is for a draw. For each variable (corresponding to the columns),
   the total count of the number of times
   that variable is used in a tree decision rule (over all trees) is given.}

   \item{sigest}{
   The rough error standard deviation (\eqn{\sigma}{sigma}) used in the prior.
   }
}
%% \references{
%% Chipman, H., George, E., and McCulloch R. (2010)
%%    Bayesian Additive Regression Trees.
%%    \emph{The Annals of Applied Statistics}, \bold{4,1}, 266-298 <doi:10.1214/09-AOAS285>.

%% Chipman, H., George, E., and McCulloch R. (2006)
%%    Bayesian Ensemble Learning.
%%    Advances in Neural Information Processing Systems 19,
%%    Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.

%% Friedman, J.H. (1991)
%%    Multivariate adaptive regression splines.
%%    \emph{The Annals of Statistics}, \bold{19}, 1--67.

%% Linero, A.R. (2018)
%%   Bayesian regression trees for high dimensional prediction and variable
%%   selection. \emph{JASA}, \bold{113}, 626--36.
%% }
%% \author{
%% Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
%% Rodney Sparapani: \email{rsparapa@mcw.edu}.
%% }
\seealso{
\code{\link{pbart}}
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
Ey = f(x)
y=Ey+sigma*rnorm(n)
lmFit = lm(y~.,data.frame(x,y)) #compare lm fit to BART later

##test BART with token run to ensure installation works
set.seed(99)
bartFit = wbart(x,y,nskip=5,ndpost=5)

\dontrun{
##run BART
set.seed(99)
bartFit = wbart(x,y)

##compare BART fit to linear matter and truth = Ey
fitmat = cbind(y,Ey,lmFit$fitted,bartFit$yhat.train.mean)
colnames(fitmat) = c('y','Ey','lm','bart')
print(cor(fitmat))
}
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
