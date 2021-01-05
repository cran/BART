\name{surv.bart}
\title{Survival analysis with BART}
\alias{surv.bart}
\alias{mc.surv.bart}
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

surv.bart( x.train=matrix(0,0,0),
           y.train=NULL, times=NULL, delta=NULL,
           x.test=matrix(0,0,0),
           K=NULL, events=NULL, ztimes=NULL, zdelta=NULL,
           sparse=FALSE, theta=0, omega=1,
           a=0.5, b=1, augment=FALSE, rho=NULL,
           xinfo=matrix(0,0,0), usequants=FALSE,
           %cont=FALSE,
           rm.const=TRUE, type='pbart',
           ntype=as.integer(
               factor(type, levels=c('wbart', 'pbart', 'lbart'))),
           k=2, power=2, base=.95,
           offset=NULL, tau.num=c(NA, 3, 6)[ntype], %binaryOffset=NULL, 
           ntree=50, numcut=100, ndpost=1000, nskip=250,
           keepevery = 10L,
           %nkeeptrain=ndpost,
           %nkeeptest=ndpost,
           %nkeeptestmean=ndpost,
           %nkeeptreedraws=ndpost,
           printevery=100L, %treesaslists=FALSE,
           %keeptrainfits=TRUE,
           id=NULL,    ## surv.bart only
           seed=99,    ## mc.surv.bart only
           mc.cores=2, ## mc.surv.bart only
           nice=19L    ## mc.surv.bart only
         )

mc.surv.bart( x.train=matrix(0,0,0),
              y.train=NULL, times=NULL, delta=NULL, 
              x.test=matrix(0,0,0), 
              K=NULL, events=NULL, ztimes=NULL, zdelta=NULL,
              sparse=FALSE, theta=0, omega=1,
              a=0.5, b=1, augment=FALSE, rho=NULL,
              xinfo=matrix(0,0,0), usequants=FALSE,
              %cont=FALSE,
              rm.const=TRUE, type='pbart',
              ntype=as.integer(
                  factor(type, levels=c('wbart', 'pbart', 'lbart'))),
              k=2, power=2, base=.95,
              offset=NULL, tau.num=c(NA, 3, 6)[ntype], %binaryOffset=NULL,
              ntree=50, numcut=100, ndpost=1000, nskip=250,
              keepevery = 10L,
              %nkeeptrain=ndpost,
              %nkeeptest=ndpost,
              %nkeeptestmean=ndpost,
              %nkeeptreedraws=ndpost,
              printevery=100L, %treesaslists=FALSE,
              %keeptrainfits=TRUE,
              id=NULL,    ## surv.bart only
              seed=99,    ## mc.surv.bart only
              mc.cores=2, ## mc.surv.bart only
              nice=19L    ## mc.surv.bart only
            )
}

\arguments{

  \item{x.train}{ Explanatory variables for training (in sample)
data.\cr Must be a matrix with (as usual) rows corresponding to
observations and columns to variables.\cr \code{surv.bart} will generate
draws of \eqn{f(t, x)} for each \eqn{x} which is a row of x.train (note
that the definition of \code{x.train} is dependent on whether
\code{y.train} has been specified; see below).  }

   \item{y.train}{
   Binary response dependent variable for training (in sample) data.\cr
   If \code{y.train} is \code{NULL}, then \code{y.train} (\code{x.train} and
   \code{x.test}, if specified) are generated
   by a call to \code{surv.pre.bart} (which require that \code{times} and \code{delta} be
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
   \code{surv.bart} will generate draws of \eqn{f(t, x)} for each \eqn{x} which is a row of x.test.
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

   \item{offset}{  With binary
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

   \item{ndpost}{
   The number of posterior draws returned.
   }

   \item{nskip}{
   Number of MCMC iterations to be treated as burn in.
   }

   \item{printevery}{
   As the MCMC runs, a message is printed every printevery draws.
   }

   \item{keepevery}{
   Every keepevery draw is kept to be returned to the user.\cr
   A \dQuote{draw} will consist of values \eqn{f^*(t, x)}{f*(t, x)}
   at \eqn{x} = rows from the train(optionally) and test data, where \eqn{f^*}{f*} denotes
   the current draw of \eqn{f}.
   }

   %% \item{keeptrainfits}{
   %% If true the draws of \eqn{f(t, x)} for \eqn{x} = rows of x.train are returned.
   %% }

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

   %% \item{treesaslists}{
   %%   The default behavior is to return the trees as a text string.  If you
   %%   specify \code{treesaslists=TRUE}, then the trees will be returned
   %%   as a list as well.
   %% }

   \item{id}{
     \code{surv.bart} only: unique identifier added to returned list.
   }

   \item{seed}{
     \code{mc.surv.bart} only: seed required for reproducible MCMC.
   }

   \item{mc.cores}{
     \code{mc.surv.bart} only: number of cores to employ in parallel.
   }

   \item{nice}{
     \code{mc.surv.bart} only: set the job niceness.  The default
     niceness is 19: niceness goes from 0 (highest) to 19 (lowest).
   }
}

\value{
  \code{surv.bart} returns an object of type \code{survbart} which is
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

   \item{yhat.test}{Same as yhat.train but now the x's are the rows of the test data.}
   \item{surv.test}{The survival function, \eqn{S(t|x)}, where x's are the rows of the test data.}
   \item{yhat.train.mean}{train data fits = mean of yhat.train columns.}
   \item{yhat.test.mean}{test data fits = mean of yhat.test columns.}
   \item{surv.test.mean}{mean of surv.test columns.}
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

%% Sparapani, R., Logan, B., McCulloch, R., and Laud, P. (2016)
%%   Nonparametric survival analysis using Bayesian Additive Regression Trees (BART).
%%   \emph{Statistics in Medicine}, \bold{16}, 2741--53 <doi:10.1002/sim.6893>.
 
%% }
%% \author{
%% Rodney Sparapani: \email{rsparapa@mcw.edu}
%% }
\seealso{
\code{\link{surv.pre.bart}}
}
\examples{

## load survival package for the advanced lung cancer example
data(lung)

N <- length(lung$status)

table(lung$ph.karno, lung$pat.karno)

## if physician's KPS unavailable, then use the patient's
h <- which(is.na(lung$ph.karno))
lung$ph.karno[h] <- lung$pat.karno[h]

times <- lung$time
delta <- lung$status-1 ##lung$status: 1=censored, 2=dead
##delta: 0=censored, 1=dead

## this study reports time in days rather than weeks or months
## coarsening from days to weeks or months will reduce the computational burden
##times <- ceiling(times/30)
times <- ceiling(times/7)  ## weeks

table(times)
table(delta)

## matrix of observed covariates
x.train <- cbind(lung$sex, lung$age, lung$ph.karno)

## lung$sex:        Male=1 Female=2
## lung$age:        Age in years
## lung$ph.karno:   Karnofsky performance score (dead=0:normal=100:by=10)
##                  rated by physician

dimnames(x.train)[[2]] <- c('M(1):F(2)', 'age(39:82)', 'ph.karno(50:100:10)')

table(x.train[ , 1])
summary(x.train[ , 2])
table(x.train[ , 3])

##test BART with token run to ensure installation works
set.seed(99)
post <- surv.bart(x.train=x.train, times=times, delta=delta,
                  nskip=1, ndpost=1, keepevery=1)

\dontrun{
## run one long MCMC chain in one process
## set.seed(99)
## post <- surv.bart(x.train=x.train, times=times, delta=delta, x.test=x.test)

## in the interest of time, consider speeding it up by parallel processing
## run "mc.cores" number of shorter MCMC chains in parallel processes
post <- mc.surv.bart(x.train=x.train, times=times, delta=delta,
                     mc.cores=8, seed=99)

pre <- surv.pre.bart(times=times, delta=delta, x.train=x.train,
                     x.test=x.train)

K <- pre$K
M <- nrow(post$yhat.train)

pre$tx.test <- rbind(pre$tx.test, pre$tx.test)
pre$tx.test[ , 2] <- c(rep(1, N*K), rep(2, N*K))
## sex pushed to col 2, since time is always in col 1

pred <- predict(post, newdata=pre$tx.test, mc.cores=8)

pd <- matrix(nrow=M, ncol=2*K)

for(j in 1:K) {
    h <- seq(j, N*K, by=K)
    pd[ , j] <- apply(pred$surv.test[ , h], 1, mean)
    pd[ , j+K] <- apply(pred$surv.test[ , h+N*K], 1, mean)
}

pd.mu  <- apply(pd, 2, mean)
pd.025 <- apply(pd, 2, quantile, probs=0.025)
pd.975 <- apply(pd, 2, quantile, probs=0.975)

males <- 1:K
females <- males+K

plot(c(0, pre$times), c(1, pd.mu[males]), type='s', col='blue',
     ylim=0:1, ylab='S(t, x)', xlab='t (weeks)',
     main=paste('Advanced Lung Cancer ex. (BART::lung)',
                "Friedman's partial dependence function",
                'Male (blue) vs. Female (red)', sep='\n'))
lines(c(0, pre$times), c(1, pd.025[males]), col='blue', type='s', lty=2)
lines(c(0, pre$times), c(1, pd.975[males]), col='blue', type='s', lty=2)
lines(c(0, pre$times), c(1, pd.mu[females]), col='red', type='s')
lines(c(0, pre$times), c(1, pd.025[females]), col='red', type='s', lty=2)
lines(c(0, pre$times), c(1, pd.975[females]), col='red', type='s', lty=2)


}
}
\keyword{nonparametric survival model}
\keyword{nonproportional hazards}
