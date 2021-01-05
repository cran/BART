\name{rs.pbart}
\title{BART for dichotomous outcomes with parallel computation and
  stratified random sampling}
\alias{rs.pbart}
\description{
BART is a Bayesian \dQuote{sum-of-trees} model.\cr
For numeric response \eqn{y}, we have
\eqn{y = f(x) + \epsilon}{y = f(x) + e},
where \eqn{\epsilon \sim N(0,\sigma^2)}{e ~ N(0,sigma^2)}.\cr
For a binary response \eqn{y}, \eqn{P(Y=1 | x) = F(f(x))}, where \eqn{F}
denotes the standard normal cdf (probit link).

In both cases, \eqn{f} is the sum of many tree models.
The goal is to have very flexible inference for the uknown
function \eqn{f}.

In the spirit of \dQuote{ensemble models},
each tree is constrained by a prior to be a weak learner
so that it contributes a
small amount to the overall fit.
}
\usage{
rs.pbart(
   x.train, y.train, x.test=matrix(0.0,0,0),
   C=floor(length(y.train)/2000),
   k=2.0, power=2.0, base=.95,
   binaryOffset=0,
   ntree=50L, numcut=100L,
   ndpost=1000L, nskip=100L,
   keepevery=1L, printevery=100,
   keeptrainfits=FALSE, transposed=FALSE,
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
   \code{pbart} will generate draws of \eqn{f(x)} for each \eqn{x}
   which is a row of x.train.
   }

   \item{y.train}{
   Dependent variable for training (in sample) data.\cr
   If y is numeric a continous response model is fit (normal errors).\cr
   If y is a factor (or just has values 0 and 1) then a binary response model
   with a probit link is fit.
   }

   \item{x.test}{
   Explanatory variables for test (out of sample) data.\cr
   Should have same structure as x.train.\cr
   \code{pbart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of x.test.
   }

   \item{C}{
   The number of shards to break the data into and analyze separately.
 }
 
   \item{k}{
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
   The idea is that \eqn{f} is shrunk towards 0, so the offset allows you to shrink towards
   a probability other than .5.
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
   When running \code{pbart} in parallel, it is more memory-efficient
   to transpose \code{x.train} and \code{x.test}, if any, prior to
   calling \code{mc.pbart}.
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
   \code{rs.pbart} returns an object of type \code{pbart} which is
   essentially a list. % assigned class \sQuote{bart}.

   \item{yhat.shard}{
   Estimates generated from the individual shards rather than from the
   whole.  This object is only useful for assessing convergence.

   A matrix with ndpost rows and nrow(x.train) columns.
   Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{f}
   and each column corresponds to a row of x.train.
   The \eqn{(i,j)} value is \eqn{f^*(x)}{f*(x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
   and the \eqn{j^{th}}{j\^th} row of x.train.\cr
   Burn-in is dropped.
 }
 
 \item{yhat.train}{
   Estimates generated from the whole if \code{keeptrainfits=TRUE}.
   
   A matrix with ndpost rows and nrow(x.train) columns.
   Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{f}
   and each column corresponds to a row of x.train.
   The \eqn{(i,j)} value is \eqn{f^*(x)}{f*(x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
   and the \eqn{j^{th}}{j\^th} row of x.train.\cr
   Burn-in is dropped.
   }

   \item{yhat.test}{
   Estimates generated from the whole if \code{x.test} is provided.

   Same as yhat.train but now the x's are the rows of the test data.}
   %\item{yhat.train.mean}{train data fits = mean of yhat.train columns.}
   %\item{yhat.test.mean}{test data fits = mean of yhat.test columns.}
   \item{varcount}{a matrix with ndpost rows and nrow(x.train) columns.
   Each row is for a draw. For each variable (corresponding to the columns),
   the total count of the number of times
   that variable is used in a tree decision rule (over all trees) is given.}

   In addition the list
   has a binaryOffset component giving the value used.

   Note that in the binary \eqn{y}, case yhat.train and yhat.test are
   \eqn{f(x)} + binaryOffset.  If you want draws of the probability
   \eqn{P(Y=1 | x)} you need to apply the normal cdf (\code{pnorm})
   to these values.
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
\code{\link{mc.pbart}}
}
\examples{

##simulate from Friedman's five-dimensional test function
##Friedman JH. Multivariate adaptive regression splines
##(with discussion and a rejoinder by the author).
##Annals of Statistics 1991; 19:1-67.

f = function(x) #only the first 5 matter
    sin(pi*x[ , 1]*x[ , 2]) + 2*(x[ , 3]-.5)^2+x[ , 4]+0.5*x[ , 5]-1.5

sigma = 1.0  #y = f(x) + sigma*z where z~N(0, 1)
k = 50       #number of covariates
thin = 25
ndpost = 2500
nskip = 100
C = 10
m = 10
n = 10000

set.seed(12)
x.train=matrix(runif(n*k), n, k)
Ey.train = f(x.train)
y.train=(Ey.train+sigma*rnorm(n)>0)*1
table(y.train)/n

x <- x.train
x4 <- seq(0, 1, length.out=m)

for(i in 1:m) {
    x[ , 4] <- x4[i]

    if(i==1) x.test <- x
    else x.test <- rbind(x.test, x)
}

## parallel::mcparallel/mccollect do not exist on windows
if(.Platform$OS.type=='unix') {
##test BART with token run to ensure installation works
    post = rs.pbart(x.train, y.train, 
                C=C, mc.cores=4, keepevery=1,
                seed=99, ndpost=1, nskip=1)
}

\dontrun{
post = rs.pbart(x.train, y.train, x.test=x.test,
                C=C, mc.cores=8, keepevery=thin,
                seed=99, ndpost=ndpost, nskip=nskip)
str(post)

par(mfrow=c(2, 2))

M <- nrow(post$yhat.test)
pred <- matrix(nrow=M, ncol=10)

for(i in 1:m) {
    h <- (i-1)*n+1:n
    pred[ , i] <- apply(pnorm(post$yhat.test[ , h]), 1, mean)
}

pred <- apply(pred, 2, mean)

plot(x4, qnorm(pred), xlab=expression(x[4]),
     ylab='partial dependence function', type='l')

i <- floor(seq(1, n, length.out=10))
j <- seq(-0.5, 0.4, length.out=10)
for(h in 1:10) {
    auto.corr <- acf(post$yhat.shard[ , i[h]], plot=FALSE)
    if(h==1) {
        max.lag <- max(auto.corr$lag[ , 1, 1])
        plot(1:max.lag+j[h], auto.corr$acf[1+(1:max.lag), 1, 1],
             type='h', xlim=c(0, max.lag+1), ylim=c(-1, 1),
             ylab='auto-correlation', xlab='lag')
    }
    else 
        lines(1:max.lag+j[h], auto.corr$acf[1+(1:max.lag), 1, 1],
              type='h', col=h)
}

for(j in 1:10) {
    if(j==1)
        plot(pnorm(post$yhat.shard[ , i[j]]),
             type='l', ylim=c(0, 1),
             sub=paste0('N:', n, ', k:', k),
             ylab=expression(Phi(f(x))), xlab='m')
    else
        lines(pnorm(post$yhat.shard[ , i[j]]),
              type='l', col=j)
}

geweke <- gewekediag(post$yhat.shard)

j <- -10^(log10(n)-1)
plot(geweke$z, pch='.', cex=2, ylab='z', xlab='i',
     sub=paste0('N:', n, ', k:', k),
     xlim=c(j, n), ylim=c(-5, 5))
lines(1:n, rep(-1.96, n), type='l', col=6)
lines(1:n, rep(+1.96, n), type='l', col=6)
lines(1:n, rep(-2.576, n), type='l', col=5)
lines(1:n, rep(+2.576, n), type='l', col=5)
lines(1:n, rep(-3.291, n), type='l', col=4)
lines(1:n, rep(+3.291, n), type='l', col=4)
lines(1:n, rep(-3.891, n), type='l', col=3)
lines(1:n, rep(+3.891, n), type='l', col=3)
lines(1:n, rep(-4.417, n), type='l', col=2)
lines(1:n, rep(+4.417, n), type='l', col=2)
text(c(1, 1), c(-1.96, 1.96), pos=2, cex=0.6, labels='0.95')
text(c(1, 1), c(-2.576, 2.576), pos=2, cex=0.6, labels='0.99')
text(c(1, 1), c(-3.291, 3.291), pos=2, cex=0.6, labels='0.999')
text(c(1, 1), c(-3.891, 3.891), pos=2, cex=0.6, labels='0.9999')
text(c(1, 1), c(-4.417, 4.417), pos=2, cex=0.6, labels='0.99999')

par(mfrow=c(1, 1))

##dev.copy2pdf(file='geweke.rs.pbart.pdf')
}
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
