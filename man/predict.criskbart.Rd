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
%%         \emph{The Annals of Statistics}, \bold{19}, 1--67.
%% }
%% \author{
%% Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
%% Rodney Sparapani: \email{rsparapa@mcw.edu}.
%% }
\seealso{
\code{\link{crisk.bart}}, \code{\link{mc.crisk.bart}}, \code{\link{mc.crisk.pwbart}}, \code{\link{mc.cores.openmp}}
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
        post <- mc.crisk.bart(x.train=x.train, times=times, delta=delta,
                               seed=99, mc.cores=2, nskip=5, ndpost=5,
                               keepevery=1)

        pre <- surv.pre.bart(x.train=x.train, x.test=x.test,
                             times=times, delta=delta)

        K <- post$K

        pred <- mc.crisk.pwbart(pre$tx.test, pre$tx.test,
                                post$treedraws, post$treedraws2,
                                post$binaryOffset, post$binaryOffset2)
}

\dontrun{

## run one long MCMC chain in one process
## set.seed(99)
## post <- crisk.bart(x.train=x.train, times=times, delta=delta, x.test=x.test)

## in the interest of time, consider speeding it up by parallel processing
## run "mc.cores" number of shorter MCMC chains in parallel processes
post <- mc.crisk.bart(x.train=x.train,
                       times=times, delta=delta,
                       x.test=x.test, seed=99, mc.cores=8)

## check <- mc.crisk.pwbart(post$tx.test, post$tx.test,
##                           post$treedraws, post$treedraws2,
##                           post$binaryOffset,
##                           post$binaryOffset2, mc.cores=8)
check <- predict(post, newdata=post$tx.test, newdata2=post$tx.test2,
                 mc.cores=8)

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
