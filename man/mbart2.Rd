\name{mbart2}

\title{Multinomial BART for categorical outcomes with more categories}
\alias{mbart2}
\alias{mc.mbart2}

\description{
  BART is a Bayesian \dQuote{sum-of-trees} model.\cr
  For numeric response \eqn{y}, we have \eqn{y = f(x) +\epsilon}
  {y = f(x) + e}, where \eqn{\epsilon \sim N(0, 1)}{e ~ N(0, 1)}.\cr
For a multinomial response \eqn{y}, \eqn{P(Y=y | x) = F(f(x))},
where \eqn{F} denotes the standard Normal CDF (probit link) or the
standard Logistic CDF (logit link).

In both cases, \eqn{f} is the sum of many tree models.
The goal is to have very flexible inference for the uknown
function \eqn{f}.

In the spirit of \dQuote{ensemble models},
each tree is constrained by a prior to be a weak learner
so that it contributes a
small amount to the overall fit.
}
\usage{
mbart2(
      x.train, y.train,
      x.test=matrix(0,0,0), type='lbart',
      ntype=as.integer(
          factor(type,
                 levels=c('wbart', 'pbart', 'lbart'))),
      sparse=FALSE, theta=0, omega=1,
      a=0.5, b=1, augment=FALSE, rho=NULL,
      xinfo=matrix(0,0,0), usequants=FALSE,
      rm.const=TRUE,
      k=2, power=2, base=0.95,
      tau.num=c(NA, 3, 6)[ntype],
      offset=NULL, 
      ntree=c(200L, 50L, 50L)[ntype], numcut=100L,
      ndpost=1000L, nskip=100L,
      keepevery=c(1L, 10L, 10L)[ntype],
      printevery=100L, transposed=FALSE,
      hostname=FALSE,
      mc.cores = 2L, ## mc.bart only
      nice = 19L,    ## mc.bart only
      seed = 99L     ## mc.bart only
     )

mc.mbart2(
         x.train, y.train,
         x.test=matrix(0,0,0), type='lbart',
         ntype=as.integer(
             factor(type,
                    levels=c('wbart', 'pbart', 'lbart'))),
         sparse=FALSE, theta=0, omega=1,
         a=0.5, b=1, augment=FALSE, rho=NULL,
         xinfo=matrix(0,0,0), usequants=FALSE,
         rm.const=TRUE,
         k=2, power=2, base=0.95,
         tau.num=c(NA, 3, 6)[ntype],
         offset=NULL, 
         ntree=c(200L, 50L, 50L)[ntype], numcut=100L,
         ndpost=1000L, nskip=100L,
         keepevery=c(1L, 10L, 10L)[ntype],
         printevery=100L, transposed=FALSE,
         hostname=FALSE,
         mc.cores = 2L, ## mc.bart only
         nice = 19L,    ## mc.bart only
         seed = 99L     ## mc.bart only
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
   \code{mbart2} will generate draws of \eqn{f(x)} for each \eqn{x}
   which is a row of \code{x.train}.
   }

   \item{y.train}{
   Categorical dependent variable for training (in sample) data.\cr
   %% If y is numeric a continous response model is fit (Logistic errors).\cr
   %% If y is a factor (or just has values 0 and 1) then a binary response model
   %% with a logit link is fit.
   }

   \item{x.test}{
   Explanatory variables for test (out of sample) data.\cr
   Should have same structure as \code{x.train}.\cr
   \code{mbart2} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of \code{x.test}.
 }

 \item{type}{ You can use this argument to specify the type of fit.
   \code{'pbart'} for probit BART or \code{'lbart'} for logit BART. }

 \item{ntype}{ The integer equivalent of \code{type} where
  \code{'pbart'} is 2 and \code{'lbart'} is 3.}

   \item{sparse}{Whether to perform variable selection based on a
     sparse Dirichlet prior rather than simply uniform; see Linero
     2016.}
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
 
   \item{rm.const}{ Whether or not to remove constant variables.}

   %% \item{tau.interval}{
   %%   The width of the interval to scale the variance for the terminal
   %%   leaf values.}
   
   \item{k}{
   %% For numeric y,
   %% k is the number of prior standard deviations \eqn{E(Y|x) = f(x)} is away from +/-.5.
   %% The response (y.train) is internally scaled to range from -.5 to .5.
   For categorical \code{y.train},
   k is the number of prior standard deviations \eqn{f(x)} is away from +/-3.
   %%In both cases, the bigger k is, the more conservative the fitting will be.
   }

   \item{power}{
   Power parameter for tree prior.
   }

   \item{base}{
   Base parameter for tree prior.
   }
 
 \item{tau.num}{ The numerator in the \code{tau} definition, i.e.,
   \code{tau=tau.num/(k*sqrt(ntree))}. }

   \item{offset}{ With Multinomial
   BART, the centering is \eqn{P(yj=1 | x) = F(fj(x) + offset[j])} where
   \code{offset} defaults to \code{F^{-1}(mean(y.train))}.  You can use
   the \code{offset} parameter to over-ride these defaults.}
   
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

   %% \item{keepevery}{
   %% Every keepevery draw is kept to be returned to the user.\cr
   %% A \dQuote{draw} will consist of values of the error standard deviation (\eqn{\sigma}{sigma})
   %% and \eqn{f^*(x)}{f*(x)}
   %% at \eqn{x} = rows from the train(optionally) and test data, where \eqn{f^*}{f*} denotes
   %% the current draw of \eqn{f}.
   %% }

   \item{transposed}{
   When running \code{mbart2} in parallel, it is more memory-efficient
   to transpose \code{x.train} and \code{x.test}, if any, prior to
   calling \code{mc.mbart2}.
   }

 \item{hostname}{
   When running on a cluster occasionally it is useful
   to track on which node each chain is running; to do so
   set this argument to \code{TRUE}.
 }
 
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
   At each MCMC interation, we produce a draw from %% the joint posterior
   %% \eqn{f | (x,y)}{f \| (x,y)} in the numeric \eqn{y} case
   %% and just
   \eqn{f} in the categorical \eqn{y} case.

   Thus, unlike a lot of other modelling methods in R, we do not produce
   a single model object from which fits and summaries may be extracted.
   The output consists of values \eqn{f^*(x)}{f*(x)}
   %(and \eqn{\sigma^*}{sigma*} in the numeric case)
   where * denotes a particular draw.
   The \eqn{x} is either a row from the training data (x.train).
}
\value{
   \code{mbart2} returns an object of type \code{mbart2} which is
   essentially a list. % assigned class \sQuote{bart}.

   \item{yhat.train}{ A matrix with \code{ndpost} rows and
   \code{nrow(x.train)*K} columns.  Each row corresponds to a draw
   \eqn{f^*}{f*} from the posterior of \eqn{f} and each column
   corresponds to an estimate for a row of \code{x.train}.  For the
   \code{i}th row of \code{x.train}, we provide the corresponding
   \code{(i-1)*K+j}th column of \code{yhat.train} where
   \code{j=1,...,K} indexes the categories.\cr
   Burn-in is dropped.  }

   \item{yhat.train.mean}{train data fits = mean of \code{yhat.train}
     columns.}
   
   \item{varcount}{a matrix with \code{ndpost} rows and
   \code{nrow(x.train)} columns.  Each row is for a draw. For each
   variable (corresponding to the columns), the total count of the
   number of times that variable is used in a tree decision rule (over
   all trees) is given.}

   In addition, the list
   has a \code{offset} vector giving the value used.

   Note that in the multinomial \eqn{y} case \code{yhat.train} is 
   \eqn{f(x) + offset[j]}.
   %% If you want draws of the probability
   %% \eqn{P(Y=y | x)} you need to apply the Normal CDF (\code{plogis})
   %% to these values.
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
%% Rodney Sparapani: \email{rsparapa@mcw.edu},\cr
%% Robert Gramacy: \email{rbg@vt.edu}.
%% }
\seealso{
\code{\link{gbart}}, \code{\link{alligator}}
}
\examples{

N=500
set.seed(12)
x1=runif(N)
x2=runif(N, max=1-x1)
x3=1-x1-x2
x.train=cbind(x1, x2, x3)
y.train=0
for(i in 1:N)
    y.train[i]=sum((1:3)*rmultinom(1, 1, x.train[i, ]))
table(y.train)/N

##test mbart2 with token run to ensure installation works
set.seed(99)
post = mbart2(x.train, y.train, nskip=1, ndpost=1)

\dontrun{
set.seed(99)
post=mbart2(x.train, y.train, x.train)
##mc.post=mbart2(x.train, y.train, x.test, mc.cores=8, seed=99)

K=3
i=seq(1, N*K, K)-1
for(j in 1:K)
    print(cor(x.train[ , j], post$prob.test.mean[i+j])^2)

}
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
