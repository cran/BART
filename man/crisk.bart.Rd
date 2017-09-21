\name{crisk.bart}
\title{BART for competing risks}
\alias{crisk.bart}
\alias{mc.crisk.bart}
\description{
Here we have implemented a simple and direct approach to utilize BART
for competing risks that is very flexible, and is akin to discrete-time
survival analysis.  Following the capabilities of BART, we allow for
maximum flexibility in modeling the dependence of competing failure times on
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

crisk.bart(x.train=matrix(0, 0L, 0L), y.train=NULL,
           x.train2=x.train, y.train2=NULL,
           times=NULL, delta=NULL,
           x.test=matrix(0.0,0,0), x.test2=x.test,
           cond=NULL, k=2.0, power=2.0, base=.95,
           binaryOffset=NULL, binaryOffset2=NULL,
           ntree=50, numcut=100, ndpost=1000, nskip=250,
           keepevery = 10L,
           nkeeptrain=ndpost,
           nkeeptest=ndpost,
           nkeeptestmean=ndpost,
           nkeeptreedraws=ndpost,
           printevery=100L, treesaslists=FALSE, keeptrainfits=TRUE,
           id=NULL,    ## crisk.bart only
           seed=99,    ## mc.crisk.bart only
           mc.cores=2, ## mc.crisk.bart only
           nice=19L    ## mc.crisk.bart only
          )

mc.crisk.bart(x.train=matrix(0, 0L, 0L), y.train=NULL,
              x.train2=x.train, y.train2=NULL,
              times=NULL, delta=NULL,
              x.test=matrix(0.0,0,0), x.test2=x.test,
              cond=NULL, k=2.0, power=2.0, base=.95,
              binaryOffset=NULL, binaryOffset2=NULL,
              ntree=50, numcut=100, ndpost=1000, nskip=250,
              keepevery = 10L,
              nkeeptrain=ndpost,
              nkeeptest=ndpost,
              nkeeptestmean=ndpost,
              nkeeptreedraws=ndpost,
              printevery=100L, treesaslists=FALSE, keeptrainfits=TRUE,
              id=NULL,    ## crisk.bart only
              seed=99,    ## mc.crisk.bart only
              mc.cores=2, ## mc.crisk.bart only
              nice=19L    ## mc.crisk.bart only
             )
}

\arguments{

  \item{x.train}{ Explanatory variables for training (in sample)
data of cause 1.\cr Must be a matrix with (as usual) rows corresponding to
observations and columns to variables.\cr \code{crisk.bart} will generate
draws of \eqn{f(t, x)} for each \eqn{x} which is a row of x.train (note
that the definition of \code{x.train} is dependent on whether
\code{y.train} has been specified; see below).  }

   \item{y.train}{
   Cause 1 binary response for training (in sample) data.\cr
   If \code{y.train} is \code{NULL}, then \code{y.train} (\code{x.train} and
   \code{x.test}, if specified) are generated
   by a call to \code{crisk.pre.bart} (which require that \code{times} and \code{delta} be
   provided: see below); otherwise, \code{y.train} (\code{x.train} and
   \code{x.test}, if specified)
   are utilized as given assuming that the data construction has already been performed.
   }

   \item{x.train2}{ Explanatory variables for training (in sample)
     data of cause 2.}

   \item{y.train2}{
   Cause 2 binary response for training (in sample) data, i.e., failure
   from any cause besides the cause of interest which is cause 1.
   Similar to \code{y.train} above.
   }

   \item{times}{
     The time of event or right-censoring.\cr
     If \code{y.train} is \code{NULL}, then \code{times} (and \code{delta})
     must be provided.
   }

   \item{delta}{
     The event indicator: 1 is a cause 1 event, 2 a cause 2 event while 0 is censored.\cr
     If \code{y.train} is \code{NULL}, then \code{delta} (and \code{times})
     must be provided.
   }

   \item{x.test}{
   Explanatory variables for test (out of sample) data of cause 1.\cr
   Must be a matrix and have the same structure as x.train.\cr
   \code{crisk.bart} will generate draws of \eqn{f(t, x)} for each \eqn{x} which is a row of x.test.
 }

    \item{x.test2}{
      Explanatory variables for test (out of sample) data of cause 2.}

   \item{cond}{A vector of indices of \code{y.train} indicating censored subjects.}
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
   Cause 1 binary offset.
   }

   \item{binaryOffset2}{
   Cause 2 binary offset.
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

   \item{keeptrainfits}{
   Whether to keep \code{yhat.train} or not.
   }

   \item{treesaslists}{
     The default behavior is to return the trees as a text string.  If you
     specify \code{treesaslists=TRUE}, then the trees will be returned
     as a list as well.
   }

   \item{id}{
     \code{crisk.bart} only: unique identifier added to returned list.
   }

   \item{seed}{
     \code{mc.crisk.bart} only: seed required for reproducible MCMC.
   }

   \item{mc.cores}{
     \code{mc.crisk.bart} only: number of cores to employ in parallel.
   }

   \item{nice}{
     \code{mc.crisk.bart} only: set the job niceness.  The default
     niceness is 19: niceness goes from 0 (highest) to 19 (lowest).
   }
}

\value{
  \code{crisk.bart} returns an object of type \code{criskbart} which is
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

   \item{yhat.test}{Same as \code{yhat.train} but now the x's are the rows of the test data.}
   \item{surv.test}{test data fits for survival probability.}
   \item{surv.test.mean}{mean of \code{surv.test} over the posterior samples.}
   \item{prob.test}{The probability of suffering cause 1 which is
   occasionally useful, e.g., in calculating the concordance.}
   \item{prob.test2}{The probability of suffering cause 2 which is
   occasionally useful, e.g., in calculating the concordance.}
   \item{cif.test}{The cumulative incidence function of cause 1,
   \eqn{F_1(t, x)}, where x's are the rows of the test data.}
   \item{cif.test2}{The cumulative incidence function of cause 2,
   \eqn{F_2(t, x)}, where x's are the rows of the test data.}
   %\item{yhat.train.mean}{train data fits = mean of yhat.train columns.}
   \item{yhat.test.mean}{test data fits = mean of yhat.test columns.}
   \item{cif.test.mean}{mean of \code{cif.test} columns for cause 1.}
   \item{cif.test2.mean}{mean of \code{cif.test2} columns for cause 2.}
   \item{varcount}{a matrix with \code{ndpost} rows and
   \code{nrow(x.train)} columns.
   Each row is for a draw. For each variable (corresponding to the columns),
   the total count of the number of times
   that variable is used for cause 1 in a tree decision rule (over all trees) is given.}
   \item{varcount2}{ For each variable the total count of the number of times
   that variable is used for cause 2 in a tree decision rule is given.}
}
\references{
  Sparapani, R., Logan, B., McCulloch, R., and Laud, P. (2016)
  Nonparametric survival analysis using Bayesian Additive Regression Trees (BART).
   \emph{Statistics in Medicine}, \bold{16:2741-53} <doi:10.1002/sim.6893>.
}
\author{
Rodney Sparapani: \email{rsparapa@mcw.edu}
}
\seealso{
\code{\link{crisk.pre.bart}}, \code{\link{predict.criskbart}},
   \code{\link{mc.crisk.pwbart}}
}
\examples{

data(transplant)

pfit <- survfit(Surv(futime, event) ~ abo, transplant)

# competing risks for type O
plot(pfit[4,], xscale=30.5, xmax=735, col=1:3, lwd=2, ylim=c(0, 1),
       xlab='t (months)', ylab='Aalen-Johansen (AJ) CI(t)')
    legend(450, .4, c("Death", "Transplant", "Withdrawal"), col=1:3, lwd=2)

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

##test BART with token run to ensure installation works
set.seed(99)
post <- crisk.bart(x.train=x.train, times=times, delta=delta, x.test=x.test,
                   nskip=5, ndpost=5, keepevery=1)

\dontrun{
## run one long MCMC chain in one process
set.seed(99)
post <- crisk.bart(x.train=x.train, times=times, delta=delta, x.test=x.test)
## in the interest of time, consider speeding it up by parallel processing
## run "mc.cores" number of shorter MCMC chains in parallel processes
## post <- mc.crisk.bart(x.train=x.train, times=times, delta=delta, x.test=x.test,
##                       seed=99, mc.cores=5)

K <- post$K

typeO.cif <- matrix(0, nrow=nrow(post$cif.test), ncol=K)

for(i in 1:N) for(j in 1:K) {
                  h <- (i-1)*K+j
                  typeO.cif[ , j] <- typeO.cif[ , j]+ post$cif.test[ , h]/N
                  }

typeO.cif.mean <- apply(typeO.cif, 2, mean)
typeO.cif.025 <- apply(typeO.cif, 2, quantile, probs=0.025)
typeO.cif.975 <- apply(typeO.cif, 2, quantile, probs=0.975)

plot(pfit[4,], xscale=30.5, xmax=735, col=1:3, lwd=2, ylim=c(0, 0.8),
       xlab='t (months)', ylab='CI(t)')
points(c(0, post$times)*30.5, c(0, typeO.cif.mean), col=4, type='s', lwd=2)
points(c(0, post$times)*30.5, c(0, typeO.cif.025), col=4, type='s', lwd=2, lty=2)
points(c(0, post$times)*30.5, c(0, typeO.cif.975), col=4, type='s', lwd=2, lty=2)
     legend(450, .4, c("Transplant(BART)", "Transplant(AJ)",
                       "Death(AJ)", "Withdrawal(AJ)"),
            col=c(4, 2, 1, 3), lwd=2)

}
}
\keyword{nonparametric survival model}
\keyword{nonproportional hazards}
