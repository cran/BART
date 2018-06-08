\name{mbart}
\title{BART for multinomial outcomes with Normal latents}
\alias{mbart}
\description{
BART is a Bayesian \dQuote{sum-of-trees} model.\cr
For numeric response \eqn{y}, we have
\eqn{y = f(x) + \epsilon}{y = f(x) + e},
where \eqn{\epsilon \sim N(0, 1)}{e ~ N(0, 1)}.\cr
For a multinomial response \eqn{y}, \eqn{P(Y=y | x) = F(f(x))}, where \eqn{F}
denotes the standard Normal CDF (probit link).

In both cases, \eqn{f} is the sum of many tree models.
The goal is to have very flexible inference for the uknown
function \eqn{f}.

In the spirit of \dQuote{ensemble models},
each tree is constrained by a prior to be a weak learner
so that it contributes a
small amount to the overall fit.
}
\usage{
mbart(
   x.train, y.train, x.test=matrix(0.0,0,0),
   sparse=FALSE, theta=0, omega=1,
   a=0.5, b=1, augment=FALSE, rho=NULL,
   xinfo=matrix(0.0,0,0), usequants=FALSE,
   cont=FALSE, rm.const=TRUE, %%tau.interval=0.95,
   k=2.0, power=2.0, base=.95, 
   binaryOffset=NULL,
   ntree=50L, numcut=100L,
   ndpost=1000L, nskip=100L,
   keepevery=1L,
   nkeeptrain=ndpost, nkeeptest=ndpost,
   %nkeeptestmean=ndpost,
   nkeeptreedraws=ndpost,
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
   \code{mbart} will generate draws of \eqn{f(x)} for each \eqn{x}
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
   \code{mbart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of \code{x.test}.
   }

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
   
   \item{cont}{ Whether or not to assume all variables are continuous.}

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

   \item{binaryOffset}{
   Used for categorical \eqn{y}.\cr
   The model is \eqn{P(Y=j | x) = F(f(x) + binaryOffset[j])}.\cr
   }
   
   %% \item{binaryOffset}{
   %% Used for binary \eqn{y}.\cr
   %% The model is \eqn{P(Y=1 | x) = F(f(x) + binaryOffset)}.\cr
   %% The idea is that \eqn{f} is shrunk towards 0, so the offset allows you to shrink towards
   %% a probability other than .5.
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

   \item{nkeeptrain}{
   Number of MCMC iterations to be returned for train data.
   }

   \item{nkeeptest}{
   Number of MCMC iterations to be returned for test data.
   }

   %% \item{nkeeptestmean}{
   %% Number of MCMC iterations to be returned for test mean.
   %% }

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

   \item{transposed}{
   When running \code{mbart} in parallel, it is more memory-efficient
   to transpose \code{x.train} and \code{x.test}, if any, prior to
   calling \code{mc.mbart}.
   }

   %% \item{treesaslists}{
   %%   The default behavior is to return the trees as a text string.  If you
   %%   specify \code{treesaslists=TRUE}, then the trees will be returned
   %%   as a list as well.
   %% }

}
\details{
   BART is an Bayesian MCMC method.
   At each MCMC interation, we produce a draw from %% the joint posterior
   %% \eqn{f | (x,y)}{f \| (x,y)} in the numeric \eqn{y} case
   %% and just
   \eqn{f} in the categorical \eqn{y} case.

   Thus, unlike a lot of other modelling methods in R, we do not produce a single model object
   from which fits and summaries may be extracted.  The output consists of values
   \eqn{f^*(x)}{f*(x)}
   %(and \eqn{\sigma^*}{sigma*} in the numeric case)
   where * denotes a particular draw.
   The \eqn{x} is either a row from the training data (x.train).
}
\value{
   \code{mbart} returns an object of type \code{mbart} which is
   essentially a list. % assigned class \sQuote{bart}.

   \item{yhat.train}{
   A matrix with ndpost rows and nrow(x.train) columns.
   Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{f}
   and each column corresponds to a row of x.train.
   The \eqn{(i,j)} value is \eqn{f^*(x)}{f*(x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
   and the \eqn{j^{th}}{j\^th} row of x.train.\cr
   Burn-in is dropped.
   }

   \item{yhat.train.mean}{train data fits = mean of yhat.train columns.}
   \item{varcount}{a matrix with ndpost rows and nrow(x.train) columns.
   Each row is for a draw. For each variable (corresponding to the columns),
   the total count of the number of times
   that variable is used in a tree decision rule (over all trees) is given.}

   In addition, the list
   has a binaryOffset vector giving the value used.

   Note that in the multinomial \eqn{y} case \code{yhat.train} is 
   \eqn{f(x) + binaryOffset[j]}.
   %% If you want draws of the probability
   %% \eqn{P(Y=y | x)} you need to apply the Normal CDF (\code{plogis})
   %% to these values.
}
\references{
Chipman, H., George, E., and McCulloch R. (2010)
   Bayesian Additive Regression Trees.
   \emph{The Annals of Applied Statistics}, \bold{4,1}, 266-298 <doi: 10.1214/09-AOAS285>.

%% Chipman, H., George, E., and McCulloch R. (2006)
%%    Bayesian Ensemble Learning.
%%    Advances in Neural Information Processing Systems 19,
%%    Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.

Friedman, J.H. (1991)
   Multivariate adaptive regression splines.
        \emph{The Annals of Statistics}, \bold{19}, 1--67.
}
\author{
Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
Rodney Sparapani: \email{rsparapa@mcw.edu},\cr
Robert Gramacy: \email{rbg@vt.edu}.
}
\seealso{
\code{\link{pbart}}
}
\examples{

##simulate from Friedman's five-dimensional test function
##Friedman JH. Multivariate adaptive regression splines
##(with discussion and a rejoinder by the author).
##Annals of Statistics 1991; 19:1-67.

f = function(x) #only the first 5 matter
    sin(pi*x[ , 1]*x[ , 2]) + 2*(x[ , 3]-.5)^2+x[ , 4]+0.5*x[ , 5]-1.5

sigma = 1.0  # y~N(f(x), sigma)
P = 5        # number of covariates
N = 2500
C = 3        # number of categories
M = 5^P

x <- seq(0, 1, length.out=5)
x.test <- matrix(nrow=M, ncol=P)
dimnames(x.test)[[2]] <- paste0('x', 1:5)
h <- 1
for(x5 in x)
    for(x4 in x)
        for(x3 in x)
            for(x2 in x)
                for(x1 in x) {
                    x.test[h, ] <- c(x1, x2, x3, x4, x5)
                    h <- h+1
                }
Ey.test = f(x.test)

set.seed(12)
x.train=matrix(runif(N*P), N, P)
dimnames(x.train)[[2]] <- paste0('x', 1:5)
Ey.train = f(x.train)
y.train=rnorm(N, Ey.train, sigma)

i <- y.train > -1
j <- y.train > 1

y.train[!i] <- 1
y.train[i] <- 2
y.train[j] <- 3

table(y.train)
##test BART with token run to ensure installation works
## set.seed(99)
## post = mbart(x.train, y.train, x.test,
##              nskip=1, ndpost=1, keepevery=1)

\dontrun{ 
## set.seed(99)
## post = mbart(x.train, y.train, x.test)

post = mc.mbart(x.train, y.train, x.test, mc.cores=8, seed=99)

h <- seq(1, C*N, by=C)

print(cor(post$prob.train.mean[h], pnorm(-1, Ey.train, sigma))^2)
print(cor(post$prob.train.mean[h+1], pnorm(1, Ey.train, sigma)-pnorm(-1, Ey.train, sigma))^2)
print(cor(post$prob.train.mean[h+2], pnorm(1, Ey.train, sigma, 0))^2)

plot(pnorm(-1, Ey.train, sigma), post$prob.train.mean[h], pch='.',
     xlim=0:1, ylim=0:1, xlab='Known P(y=1)', ylab='Est. P(y=1)')
abline(0, 1)

plot(pnorm(1, Ey.train, sigma)-pnorm(-1, Ey.train, sigma), post$prob.train.mean[h+1], pch='.',
     xlim=0:1, ylim=0:1, xlab='Known P(y=2)', ylab='Est. P(y=2)')
abline(0, 1)

plot(pnorm(1, Ey.train, sigma, 0), post$prob.train.mean[h+2], pch='.',
     xlim=0:1, ylim=0:1, xlab='Known P(y=3)', ylab='Est. P(y=3)')
abline(0, 1)

h <- seq(1, C*M, by=C)

print(cor(post$prob.test.mean[h], pnorm(-1, Ey.test, sigma))^2)
print(cor(post$prob.test.mean[h+1], pnorm(1, Ey.test, sigma)-pnorm(-1, Ey.test, sigma))^2)
print(cor(post$prob.test.mean[h+2], pnorm(1, Ey.test, sigma, 0))^2)

plot(pnorm(-1, Ey.test, sigma), post$prob.test.mean[h], pch='.',
     xlim=0:1, ylim=0:1, xlab='Known P(y=1)', ylab='Est. P(y=1)')
abline(0, 1)

plot(pnorm(1, Ey.test, sigma)-pnorm(-1, Ey.test, sigma), post$prob.test.mean[h+1], pch='.',
     xlim=0:1, ylim=0:1, xlab='Known P(y=2)', ylab='Est. P(y=2)')
abline(0, 1)

plot(pnorm(1, Ey.test, sigma, 0), post$prob.test.mean[h+2], pch='.',
     xlim=0:1, ylim=0:1, xlab='Known P(y=3)', ylab='Est. P(y=3)')
abline(0, 1)

}
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
