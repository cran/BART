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

recur.bart(x.train=matrix(0.0, 0L, 0L),
           y.train=NULL, times=NULL, delta=NULL,
           x.test=matrix(0.0,0,0), x.test.nogrid=FALSE,
           k=2.0, power=2.0, base=0.95, binaryOffset=NULL,
           ntree=50, numcut = 100L, ndpost=10000, nskip=250,
           keepevery=10, nkeeptrain = ndpost\%/\%keepevery, 
           nkeeptest=ndpost\%/\%keepevery, nkeeptestmean=ndpost\%/\%keepevery, 
           nkeeptreedraws = ndpost\%/\%keepevery,
           printevery = 100L, treesaslists = FALSE, keeptrainfits = TRUE,
           seed=99,    ## mc.recur.bart only
           mc.cores=2, ## mc.recur.bart only
           nice=19L    ## mc.recur.bart only
         )

mc.recur.bart(x.train=matrix(0.0, 0L, 0L),
              y.train=NULL, times=NULL, delta=NULL,
              x.test=matrix(0.0,0,0), x.test.nogrid=FALSE,
              k=2.0, power=2.0, base=0.95, binaryOffset=NULL,
              ntree=50, numcut = 100L, ndpost=10000, nskip=250,
              keepevery=10, nkeeptrain = ndpost\%/\%keepevery, 
              nkeeptest=ndpost\%/\%keepevery, nkeeptestmean=ndpost\%/\%keepevery, 
              nkeeptreedraws = ndpost\%/\%keepevery,
              printevery = 100L, treesaslists = FALSE, keeptrainfits = TRUE, 
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

   \item{binaryOffset}{
     The model is \eqn{P(Y=1 | t, x) = F(f(t, x) + \mu_0)}{P(Y=1 | t, x) = F(f(t, x) + mu0)}
     where \eqn{\mu_0}{mu0} is specified by \code{binaryOffset}.\cr
   The idea is that \eqn{f} is shrunk towards 0, so the offset allows you to shrink towards
   a probability other than .5.\cr
   If \code{binaryOffset=NULL} when \code{times} and \code{delta} were provided, then
   an exponential distribution offset is assumed independent of the covariates, i.e.
   \code{binaryOffset=qnorm(1-exp(-mean.diff*sum(delta)/sum(times)))} where
   \code{mean.diff} is the mean of the differences of the distinct
   ordered adjacent \code{times},
   i.e. \eqn{mean(t_{(1)}-t_{(0)}, ..., t_{(K)}-t_{(K-1)})}{mean(t(1)-t(0), ..., t(K)-t(K-1))}.\cr
   If \code{binaryOffset=NULL} when \code{times} and \code{delta} were not provided, then
   \code{binaryOffset=0}.
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
   The number of posterior draws after burn in, ndpost/keepevery will actually be returned.
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
 
   \item{treesaslists}{
     The default behavior is to return the trees as a text string.  If you
     specify \code{treesaslists=TRUE}, then the trees will be returned
     as a list as well.
   }

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
  \code{recur.bart} returns a list.  Besides the items listed
   below, the list has a \code{binaryOffset} component giving the value
   used, a \code{times} component giving the unique times, \code{K}
   which is the number of unique times, \code{tx.train} and
  \code{tx.test}, if any.

   \item{yhat.train}{
   A matrix with (ndpost/keepevery) rows and nrow(x.train) columns.
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
   \item{varcount}{a matrix with (ndpost/keepevery) rows and nrow(x.train) columns.
   Each row is for a draw. For each variable (corresponding to the columns),
   the total count of the number of times
   that variable is used in a tree decision rule (over all trees) is given.}

   Note that yhat.train and yhat.test are
   \eqn{f(t, x)} + \code{binaryOffset}.  If you want draws of the probability
   \eqn{P(Y=1 | t, x)} you need to apply the normal cdf (\code{pnorm})
   to these values.
}
\references{
Sparapani, Rein, Tarima, Jackson, Meurer (2017).  Nonparametric recurrent
events analysis with BART and an application to the hospital admissions
of patients with diabetes.  MCW Biostatistics Technical Report 64
\url{http://www.mcw.edu/FileLibrary/Groups/Biostatistics/TechReports/TechReports5175/tr064.pdf}

Xu, Daniels and Winterstein.  Sequential BART for imputation of missing
covariates.  Biostatistics 2016 doi: 10.1093/biostatistics/kxw009
with supplement
\url{http://biostatistics.oxfordjournals.org/content/early/2016/03/15/biostatistics.kxw009/suppl/DC1}
 }
\author{
Rodney Sparapani: \email{rsparapa@mcw.edu}
}
\seealso{
\code{\link{recur.pre.bart}}
}
\examples{
    
## load 20 percent random sample
data(xdm20.train)
data(xdm20.test)
data(ydm20.train)

##test BART with token run to ensure installation works
set.seed(99)    
post <- recur.bart(x.train=xdm20.train, y.train=ydm20.train,
                   nskip=5, ndpost=5, keepevery=1)
    
\dontrun{
set.seed(99)    
post <- recur.bart(x.train=xdm20.train, y.train=ydm20.train)
## larger data sets can take some time so
## if parallel processing is available, uncomment this block    
## post <- mc.recur.bart(x.train=xdm20.train, y.train=ydm20.train,
##                      mc.cores=5, seed=99)
    
require(rplot)
require(rpart.plot)
    
dss <- rpart(post$yhat.train.mean~xdm20.train)

rpart.plot(dss)
## notice that all splits in the tree
## (except 1 at the bottom involving a small group)
## involve ca, cci_pud, cci_pvd, ins270 and n

## compare patients treated with insulin vs not   
N.train <- 50
N.test <- 50
K <- post$K ## 798 unique time points

## only testing set, i.e., remove training set    
xdm20.test. <- xdm20.test[N.train*K+(1:(N.test*K)), ]
xdm20.test. <- rbind(xdm20.test., xdm20.test.)
xdm20.test.[ , 'ins270'] <- rep(0:1, each=N.test*K)

pred <- pwbart(xdm20.test., post$treedraws)
## larger data sets such as these can take some time
## if parallel processing is available, uncomment this block    
## pred <- mc.pwbart(x.test., post$treedraws, mc.cores=5)

pred.haz.test <- pnorm(pred)    

## create Friedman's partial dependence function for the
## intensity/hazard by time and ins270
NK.test <- N.test*K
M <- nrow(pred.haz.test) ## number of MCMC samples, typically 1000
    
RI <- matrix(0, M, K)
    
for(i in 1:N.test) 
    RI <- RI+(pred.haz.test[ , (N.test+i-1)*K+1:K]/
              pred.haz.test[ , (i-1)*K+1:K])/N.test
    
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
