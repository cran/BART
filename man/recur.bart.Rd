\name{recur.bart}
\title{BART for recurrent events}
\alias{recur.bart}
\alias{mc.recur.bart}
\description{
Here we have implemented a simple and direct approach to utilize BART in
survival analysis that is very flexible, and is akin to discrete-time
survival analysis.  Following the capabilities of BART, we allow for
maximum flexibility in modeling the dependence of survival times on
covariates.  In particular, we do not impose proportional hazards.

To elaborate, consider data in the usual form:
\eqn{(t_i, \delta_i, {x}_i)}{(t, delta, x)} where \eqn{t_i}{t} is the event time,
\eqn{\delta_i}{delta} is an indicator distinguishing events
(\eqn{\delta=1}{delta=1}) from right-censoring
(\eqn{\delta=0}{delta=0}), \eqn{{x}_i}{x} is a vector of covariates, and
\eqn{i=1, ..., N}{i=1, ..., N (i suppressed for convenience)}
indexes subjects.

We denote the \eqn{K} distinct event/censoring times by
\eqn{0<t_{(1)}<...<t_{(K)}<\infty}{0<t(1)<...< t(K)<infinity} thus
taking \eqn{t_{(j)}}{t(j)} to be the \eqn{j^{th}}{j'th} order
statistic among distinct observation times and, for convenience,
\eqn{t_{(0)}=0}{t(0)=0}. Now consider event indicators \eqn{y_{ij}}{y(j)}
for each subject \eqn{i} at each distinct time \eqn{t_{(j)}}{t(j)}
up to and including the subject's observation time
\eqn{t_i=t_{(n_i)}}{t=t(n)} with
\eqn{n_i=\sum_j I[t_{(j)}\leq t_i]}{n=sum I[t(j)<=t]}.
This means \eqn{y_{ij}=0}{y(j)=0} if \eqn{j<n_i}{j<n} and
\eqn{y_{in_i}=\delta_i}{y(n)=delta}.

We then denote by \eqn{p_{ij}}{p(j)} the probability
of an event at time \eqn{t_{(j)}}{t(j)} conditional on no previous event. We
now write the model for \eqn{y_{ij}}{y(j)} as a nonparametric probit
regression of \eqn{y_{ij}}{y(j)} on the time \eqn{t_{(j)}}{t(j)} and the covariates
\eqn{{x}_i}{x}, and then utilize BART for binary responses.  Specifically,
\eqn{ y_{ij}\ =\ \delta_i I[t_i=t_{(j)}],\ j=1, ..., n_i }{ y(j) = delta I[t=t(j)], j=1, ..., n }; we have
\eqn{p_{ij} = F(\mu_{ij}),\ \mu_{ij} = \mu_0+f(t_{(j)}, {x}_i)}{p(j) =
  F(mu(j)), mu(j) = mu0+f(t(j), x)} where \eqn{F} denotes the standard normal cdf (probit link).
As in the binary
response case, \eqn{f} is the sum of many tree models.  }

\usage{

recur.bart(x.train=matrix(0,0,0),
           y.train=NULL, times=NULL, delta=NULL,
           x.test=matrix(0,0,0), x.test.nogrid=FALSE,
           sparse=FALSE, theta=0, omega=1,
           a=0.5, b=1, augment=FALSE, rho=NULL,
           xinfo=matrix(0,0,0), usequants=FALSE,
           %cont=FALSE,
           rm.const=TRUE, type='pbart',
           ntype=as.integer(
               factor(type, levels=c('wbart', 'pbart', 'lbart'))),
           k=2, power=2, base=0.95,
           offset=NULL, tau.num=c(NA, 3, 6)[ntype], %binaryOffset=NULL, 
           ntree=50, numcut = 100L, ndpost=1000, nskip=250,
           keepevery=10, %nkeeptrain = ndpost,
           %nkeeptest=ndpost, %nkeeptestmean=ndpost,
           %nkeeptreedraws = ndpost,
           printevery = 100L, %treesaslists = FALSE,
           keeptrainfits = TRUE,
           seed=99,    ## mc.recur.bart only
           mc.cores=2, ## mc.recur.bart only
           nice=19L    ## mc.recur.bart only
         )

mc.recur.bart(x.train=matrix(0,0,0),
              y.train=NULL, times=NULL, delta=NULL,
              x.test=matrix(0,0,0), x.test.nogrid=FALSE,
              sparse=FALSE, theta=0, omega=1,
              a=0.5, b=1, augment=FALSE, rho=NULL,
              xinfo=matrix(0,0,0), usequants=FALSE,
              %cont=FALSE,
              rm.const=TRUE, type='pbart',
              ntype=as.integer(
                  factor(type, levels=c('wbart', 'pbart', 'lbart'))),
              k=2, power=2, base=0.95,
              offset=NULL, tau.num=c(NA, 3, 6)[ntype], %binaryOffset=NULL,
              ntree=50, numcut = 100L, ndpost=1000, nskip=250,
              keepevery=10, %nkeeptrain = ndpost,
              %nkeeptest=ndpost, %nkeeptestmean=ndpost,
              %nkeeptreedraws = ndpost,
              printevery = 100L, %treesaslists = FALSE, 
              keeptrainfits = TRUE,
              seed=99,    ## mc.recur.bart only
              mc.cores=2, ## mc.recur.bart only
              nice=19L    ## mc.recur.bart only
            )
}

\arguments{

  \item{x.train}{ Explanatory variables for training (in sample)
data.\cr Must be a matrix with (as usual) rows corresponding to
observations and columns to variables.\cr \code{recur.bart} will generate
draws of \eqn{f(t, x)} for each \eqn{x} which is a row of x.train (note
that the definition of \code{x.train} is dependent on whether
\code{y.train} has been specified; see below).  }

   \item{y.train}{
   Binary response dependent variable for training (in sample) data.\cr
   If \code{y.train} is \code{NULL}, then \code{y.train} (\code{x.train} and
   \code{x.test}, if specified) are generated
   by a call to \code{recur.pre.bart} (which require that \code{times} and \code{delta} be
   provided: see below); otherwise, \code{y.train} (\code{x.train} and
   \code{x.test}, if specified)
   are utilized as given assuming that the data construction has already been performed.
   }

   \item{times}{
     The time of event or right-censoring.\cr
     If \code{y.train} is \code{NULL}, then \code{times} (and \code{delta})
     must be provided.
   }

   \item{delta}{
     The event indicator: 1 is an event while 0 is censored.\cr
     If \code{y.train} is \code{NULL}, then \code{delta} (and \code{times})
     must be provided.
   }

   \item{x.test}{
   Explanatory variables for test (out of sample) data.\cr
   Must be a matrix and have the same structure as x.train.\cr
   \code{recur.bart} will generate draws of \eqn{f(t, x)} for each \eqn{x} which is a row of x.test.
   }

   \item{x.test.nogrid}{
     Occasionally, you do not need the entire time grid for
     \code{x.test}.  If so, then for performance reasons, you can set this argument to \code{TRUE}.
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
   
   %\item{cont}{ Whether or not to assume all variables are continuous.}

   \item{rm.const}{ Whether or not to remove constant variables.}

   \item{type}{ Whether to employ Albert-Chib, \code{'pbart'}, or
     Holmes-Held, \code{'lbart'}. }
 \item{ntype}{ The integer equivalent of \code{type} where
  \code{'wbart'} is 1, \code{'pbart'} is 2 and
  \code{'lbart'} is 3.}
   
   \item{k}{
   k is the number of prior standard deviations \eqn{f(t, x)} is away from +/-3.
   The bigger k is, the more conservative the fitting will be.
   }

   \item{power}{
   Power parameter for tree prior.
   }

   \item{base}{
   Base parameter for tree prior.
 }
 
  \item{offset}{ With binary
   BART, the centering is \eqn{P(Y=1 | x) = F(f(x) + offset)} where
   \code{offset} defaults to \code{F^{-1}(mean(y.train))}.  You can use
   the \code{offset} parameter to over-ride these defaults.}

   \item{tau.num}{ The numerator in the \code{tau} definition, i.e.,
   \code{tau=tau.num/(k*sqrt(ntree))}. }   

   %% \item{binaryOffset}{
   %%   The model is \eqn{P(Y=1 | t, x) = F(f(t, x) + \mu_0)}{P(Y=1 | t, x) = F(f(t, x) + mu0)}
   %%   where \eqn{\mu_0}{mu0} is specified by \code{binaryOffset}.\cr
   %% The idea is that \eqn{f} is shrunk towards 0, so the offset allows you to shrink towards
   %% a probability other than .5.\cr
   %% If \code{binaryOffset=NULL} when \code{times} and \code{delta} were provided, then
   %% an exponential distribution offset is assumed independent of the covariates, i.e.
   %% \code{binaryOffset=qnorm(1-exp(-mean.diff*sum(delta)/sum(times)))} where
   %% \code{mean.diff} is the mean of the differences of the distinct
   %% ordered adjacent \code{times},
   %% i.e. \eqn{mean(t_{(1)}-t_{(0)}, ..., t_{(K)}-t_{(K-1)})}{mean(t(1)-t(0), ..., t(K)-t(K-1))}.\cr
   %% If \code{binaryOffset=NULL} when \code{times} and \code{delta} were not provided, then
   %% \code{binaryOffset=0}.
   %% }

   %% \item{M}{
   %%   For large data sets, to efficiently sample \eqn{f} with Modified LISA MCMC
   %%   divide your data into \eqn{M} shards.
   %%   Defaults to 1 shard, i.e., plain BART.  See Entezari et al.
   %% }

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

   %% \item{nkeeptrain}{
   %% Number of MCMC iterations to be returned for train data.
   %% }

   %% \item{nkeeptest}{
   %% Number of MCMC iterations to be returned for test data.
   %% }

   %% %% \item{nkeeptestmean}{
   %% %% Number of MCMC iterations to be returned for test mean.
   %% %% }

   %% \item{nkeeptreedraws}{
   %% Number of MCMC iterations to be returned for tree draws.
   %% }

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

   \item{keeptrainfits}{
   Whether to keep \code{yhat.train} or not.
   }

   %% \item{treesaslists}{
   %%   The default behavior is to return the trees as a text string.  If you
   %%   specify \code{treesaslists=TRUE}, then the trees will be returned
   %%   as a list as well.
   %% }

   \item{seed}{
     \code{mc.recur.bart} only: seed required for reproducible MCMC.
   }

   \item{mc.cores}{
     \code{mc.recur.bart} only: number of cores to employ in parallel.
   }

   \item{nice}{
     \code{mc.recur.bart} only: set the job niceness.  The default
     niceness is 19: niceness goes from 0 (highest) to 19 (lowest).
   }
}

\value{
  \code{recur.bart} returns an object of type \code{recurbart} which is
  essentially a list.  Besides the items listed
   below, the list has a \code{binaryOffset} component giving the value
   used, a \code{times} component giving the unique times, \code{K}
   which is the number of unique times, \code{tx.train} and
  \code{tx.test}, if any.

   \item{yhat.train}{
   A matrix with ndpost rows and nrow(x.train) columns.
   Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{f}
   and each column corresponds to a row of x.train.
   The \eqn{(i,j)} value is \eqn{f^*(t, x)}{f*(t, x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
   and the \eqn{j^{th}}{j\^th} row of x.train.\cr
   Burn-in is dropped.
   }

   \item{haz.train}{The hazard function, \eqn{h(t|x)}, where x's are the
   rows of the training data.}
   \item{cum.train}{The cumulative hazard function, \eqn{h(t|x)}, where x's are the
   rows of the training data.}
   \item{yhat.test}{Same as yhat.train but now the x's are the rows of the test data.}
   \item{haz.test}{The hazard function, \eqn{h(t|x)}, where x's are the
   rows of the test data.}
   \item{cum.test}{The cumulative hazard function, \eqn{h(t|x)}, where x's are the
   rows of the test data.}
   \item{varcount}{a matrix with ndpost rows and nrow(x.train) columns.
   Each row is for a draw. For each variable (corresponding to the columns),
   the total count of the number of times
   that variable is used in a tree decision rule (over all trees) is given.}

   Note that yhat.train and yhat.test are
   \eqn{f(t, x)} + \code{binaryOffset}.  If you want draws of the probability
   \eqn{P(Y=1 | t, x)} you need to apply the normal cdf (\code{pnorm})
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
  
%% Sparapani RA, LE Rein, SS Tarima, TA Jackson, JR Meurer (2018).
%%   Non-parametric recurrent events analysis with BART and an application to
%%   the hospital admissions of patients with diabetes.
%%   \emph{Biostatistics}, doi:10.1093/biostatistics/kxy032
%% }

%% \author{
%% Rodney Sparapani: \email{rsparapa@mcw.edu}
%% }
\seealso{
\code{\link{recur.pre.bart}}, \code{\link{predict.recurbart}},
\code{\link{recur.pwbart}}, \code{\link{mc.recur.pwbart}}
}
\examples{

## load 20 percent random sample
data(xdm20.train)
data(xdm20.test)
data(ydm20.train)

##test BART with token run to ensure installation works
## with current technology even a token run will violate CRAN policy
## set.seed(99)
## post <- recur.bart(x.train=xdm20.train, y.train=ydm20.train,
##                    nskip=1, ndpost=1, keepevery=1)

\dontrun{

## set.seed(99)
## post <- recur.bart(x.train=xdm20.train, y.train=ydm20.train,
##                    keeptrainfits=TRUE)

## larger data sets can take some time so, if parallel processing
## is available, submit this statement instead
post <- mc.recur.bart(x.train=xdm20.train, y.train=ydm20.train,
                      keeptrainfits=TRUE, mc.cores=8, seed=99)

require(rpart)
require(rpart.plot)

post$yhat.train.mean <- apply(post$yhat.train, 2, mean)
dss <- rpart(post$yhat.train.mean~xdm20.train)

rpart.plot(dss)
## for the 20 percent sample, notice that the top splits
## involve cci_pvd and n
## for the full data set, notice that all splits
## involve ca, cci_pud, cci_pvd, ins270 and n
## (except one at the bottom involving a small group)

## compare patients treated with insulin (ins270=1) vs
## not treated with insulin (ins270=0)
N <- 50 ## 50 training patients and 50 validation patients
K <- post$K ## 798 unique time points
NK <- 50*K

## only testing set, i.e., remove training set
xdm20.test. <- xdm20.test[NK+1:NK, post$rm.const]
xdm20.test. <- rbind(xdm20.test., xdm20.test.)
xdm20.test.[ , 'ins270'] <- rep(0:1, each=NK)

## multiple threads will be utilized if available
pred <- predict(post, xdm20.test., mc.cores=8)

## create Friedman's partial dependence function for the
## relative intensity for ins270 by time
M <- nrow(pred$haz.test) ## number of MCMC samples
RI <- matrix(0, M, K)
for(j in 1:K) {
    h <- seq(j, NK, by=K)
    RI[ , j] <- apply(pred$haz.test[ , h+NK]/
                      pred$haz.test[ , h], 1, mean)
}

RI.lo <- apply(RI, 2, quantile, probs=0.025)
RI.mu <- apply(RI, 2, mean)
RI.hi <- apply(RI, 2, quantile, probs=0.975)

plot(post$times, RI.hi, type='l', lty=2, log='y',
     ylim=c(min(RI.lo, 1/RI.hi), max(1/RI.lo, RI.hi)),
     xlab='t', ylab='RI(t, x)',
     sub='insulin(ins270=1) vs. no insulin(ins270=0)',
     main='Relative intensity of hospital admissions for diabetics')
lines(post$times, RI.mu)
lines(post$times, RI.lo, lty=2)
lines(post$times, rep(1, K), col='darkgray')

## RI for insulin therapy seems fairly constant with time
mean(RI.mu)

}
}
\keyword{nonparametric recurrent events model}
\keyword{nonproportional hazards}
