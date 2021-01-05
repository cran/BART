\name{mc.crisk2.pwbart}
\title{Predicting new observations with a previously fitted BART model}
\alias{mc.crisk2.pwbart}
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
mc.crisk2.pwbart( x.test, x.test2,
                 treedraws, treedraws2,
                 binaryOffset=0, binaryOffset2=0,
                 mc.cores=2L, type='pbart',
                 transposed=FALSE, nice=19L
               )
}
\arguments{
   \item{x.test}{
   Matrix of covariates to predict \eqn{y} for cause 1.
   }

   \item{x.test2}{
   Matrix of covariates to predict \eqn{y} for cause 2.
   }

   \item{treedraws}{
     \code{$treedraws} for cause 1.
   }

   \item{treedraws2}{
     \code{$treedraws} for cause 2.
   }

   \item{binaryOffset}{
   Mean to add on to \eqn{y} prediction for cause 1.
   }
   \item{binaryOffset2}{
   Mean to add on to \eqn{y} prediction for cause 2.
   }

   \item{mc.cores}{
   Number of threads to utilize.
   }

   \item{type}{ Whether to employ Albert-Chib, \code{'pbart'}, or
  Holmes-Held, \code{'lbart'}. }

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
Returns an object of type \code{crisk2bart} which is essentially a list with components:

   \item{yhat.test}{
   A matrix with ndpost rows and nrow(x.test) columns.
   Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{f}
   and each column corresponds to a row of x.train.
   The \eqn{(i,j)} value is \eqn{f^*(x)}{f*(x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
   and the \eqn{j^{th}}{j\^th} row of x.train.\cr
   Burn-in is dropped.
   }

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

}
%% \references{
%%   Sparapani, R., Logan, B., McCulloch, R., and Laud, P. (2016)
%%   Nonparametric survival analysis using Bayesian Additive Regression Trees (BART).
%%    \emph{Statistics in Medicine}, \bold{16:2741-53} <doi:10.1002/sim.6893>.
%% }
%% \author{
%% Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
%% Rodney Sparapani: \email{rsparapa@mcw.edu}.
%% }
\seealso{
\code{\link{pwbart}}, \code{\link{crisk2.bart}}, \code{\link{mc.crisk2.bart}}
}
\examples{

data(transplant)

delta <- (as.numeric(transplant$event)-1)
## recode so that delta=1 is cause of interest; delta=2 otherwise
delta[delta==1] <- 4
delta[delta==2] <- 1
delta[delta>1] <- 2
table(delta, transplant$event)

times <- pmax(1, ceiling(transplant$futime/7)) ## weeks
##times <- pmax(1, ceiling(transplant$futime/30.5)) ## months
table(times)

typeO <- 1*(transplant$abo=='O')
typeA <- 1*(transplant$abo=='A')
typeB <- 1*(transplant$abo=='B')
typeAB <- 1*(transplant$abo=='AB')
table(typeA, typeO)

x.train <- cbind(typeO, typeA, typeB, typeAB)

x.test <- cbind(1, 0, 0, 0)
dimnames(x.test)[[2]] <- dimnames(x.train)[[2]]

## parallel::mcparallel/mccollect do not exist on windows
if(.Platform$OS.type=='unix') {
##test BART with token run to ensure installation works
        post <- mc.crisk2.bart(x.train=x.train, times=times, delta=delta,
                               seed=99, mc.cores=2, nskip=5, ndpost=5,
                               keepevery=1)

        pre <- surv.pre.bart(x.train=x.train, x.test=x.test,
                             times=times, delta=delta)

        K <- post$K

        pred <- mc.crisk2.pwbart(pre$tx.test, pre$tx.test,
                                post$treedraws, post$treedraws2,
                                post$binaryOffset, post$binaryOffset2)
}

\dontrun{

## run one long MCMC chain in one process
## set.seed(99)
## post <- crisk2.bart(x.train=x.train, times=times, delta=delta, x.test=x.test)

## in the interest of time, consider speeding it up by parallel processing
## run "mc.cores" number of shorter MCMC chains in parallel processes
post <- mc.crisk2.bart(x.train=x.train,
                       times=times, delta=delta,
                       x.test=x.test, seed=99, mc.cores=8)

check <- mc.crisk2.pwbart(post$tx.test, post$tx.test,
                          post$treedraws, post$treedraws2,
                          post$binaryOffset,
                          post$binaryOffset2, mc.cores=8)
## check <- predict(post, newdata=post$tx.test, newdata2=post$tx.test2,
##                  mc.cores=8)

print(c(post$surv.test.mean[1], check$surv.test.mean[1],
        post$surv.test.mean[1]-check$surv.test.mean[1]), digits=22)

print(all(round(post$surv.test.mean, digits=9)==
    round(check$surv.test.mean, digits=9)))

print(c(post$cif.test.mean[1], check$cif.test.mean[1],
        post$cif.test.mean[1]-check$cif.test.mean[1]), digits=22)

print(all(round(post$cif.test.mean, digits=9)==
    round(check$cif.test.mean, digits=9)))

print(c(post$cif.test2.mean[1], check$cif.test2.mean[1],
        post$cif.test2.mean[1]-check$cif.test2.mean[1]), digits=22)

print(all(round(post$cif.test2.mean, digits=9)==
    round(check$cif.test2.mean, digits=9)))


}
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
