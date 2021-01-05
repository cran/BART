\name{lbart}
\title{Logit BART for dichotomous outcomes with Logistic latents}
\alias{lbart}
\description{
BART is a Bayesian \dQuote{sum-of-trees} model.\cr
For numeric response \eqn{y}, we have
\eqn{y = f(x) + \epsilon}{y = f(x) + e},
where \eqn{\epsilon \sim Log(0, 1)}{e ~ Log(0, 1)}.\cr
For a binary response \eqn{y}, \eqn{P(Y=1 | x) = F(f(x))}, where \eqn{F}
denotes the standard Logistic CDF (logit link).

In both cases, \eqn{f} is the sum of many tree models.
The goal is to have very flexible inference for the uknown
function \eqn{f}.

In the spirit of \dQuote{ensemble models},
each tree is constrained by a prior to be a weak learner
so that it contributes a
small amount to the overall fit.
}
\usage{
lbart(
   x.train, y.train, x.test=matrix(0.0,0,0),
   sparse=FALSE, a=0.5, b=1, augment=FALSE, rho=NULL,
   xinfo=matrix(0.0,0,0), usequants=FALSE,
   cont=FALSE, rm.const=TRUE, tau.interval=0.95,
   k=2.0, power=2.0, base=.95, 
   binaryOffset=NULL,
   ntree=200L, numcut=100L,
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
   \code{lbart} will generate draws of \eqn{f(x)} for each \eqn{x}
   which is a row of x.train.
   }

   \item{y.train}{
   Binary dependent variable for training (in sample) data.\cr
   %% If y is numeric a continous response model is fit (Logistic errors).\cr
   %% If y is a factor (or just has values 0 and 1) then a binary response model
   %% with a logit link is fit.
   }

   \item{x.test}{
   Explanatory variables for test (out of sample) data.\cr
   Should have same structure as x.train.\cr
   \code{lbart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of x.test.
   }

   \item{sparse}{Whether to perform variable selection based on a
     sparse Dirichlet prior rather than simply uniform; see Linero 2016.}
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

   \item{tau.interval}{
     The width of the interval to scale the variance for the terminal
     leaf values.}
   
   \item{k}{
   For numeric y,
   k is the number of prior standard deviations \eqn{E(Y|x) = f(x)} is away from +/-.5.
   The response (y.train) is internally scaled to range from -.5 to .5.
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
   %% The idea is that \eqn{f} is shrunk towards 0, so the offset allows you to shrink towards
   %% a probability other than .5.
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
   When running \code{lbart} in parallel, it is more memory-efficient
   to transpose \code{x.train} and \code{x.test}, if any, prior to
   calling \code{mc.lbart}.
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
   \eqn{f | (x,y)}{f \| (x,y)} in the numeric \eqn{y} case
   and just \eqn{f} in the binary \eqn{y} case.

   Thus, unlike a lot of other modelling methods in R, we do not produce a single model object
   from which fits and summaries may be extracted.  The output consists of values
   \eqn{f^*(x)}{f*(x)}
   %(and \eqn{\sigma^*}{sigma*} in the numeric case)
   where * denotes a particular draw.
   The \eqn{x} is either a row from the training data (x.train) or the test data (x.test).
}
\value{
   \code{lbart} returns an object of type \code{lbart} which is
   essentially a list. % assigned class \sQuote{bart}.

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
   \item{varcount}{a matrix with ndpost rows and nrow(x.train) columns.
   Each row is for a draw. For each variable (corresponding to the columns),
   the total count of the number of times
   that variable is used in a tree decision rule (over all trees) is given.}
   In addition, the list
   has a \code{binaryOffset} giving the value used.

   Note that in the binary \eqn{y}, case yhat.train and yhat.test are
   \eqn{f(x) + binaryOffset}.
   If you want draws of the probability
   \eqn{P(Y=1 | x)} you need to apply the Logistic CDF (\code{plogis})
   to these values.
}
%% \references{
%% Chipman, H., George, E., and McCulloch R. (2010)
%%    Bayesian Additive Regression Trees.
%%    \emph{The Annals of Applied Statistics}, \bold{4,1}, 266-298 <doi: 10.1214/09-AOAS285>.

%% Chipman, H., George, E., and McCulloch R. (2006)
%%    Bayesian Ensemble Learning.
%%    Advances in Neural Information Processing Systems 19,
%%    Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.

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
%% Robert Gramacy: \email{rbg@vt.edu},\cr
%% Rodney Sparapani: \email{rsparapa@mcw.edu}.
%% }
\seealso{
\code{\link{wbart}}
}
\examples{

data(ACTG175)

## exclude those who do not have CD4 count at 96 weeks
ex <- is.na(ACTG175$cd496)
table(ex)

## inclusion criteria are CD4 counts between 200 and 500
ACTG175$cd40 <- min(500, max(250, ACTG175$cd40))

## calculate relative CD4 decline
y <- ((ACTG175$cd496-ACTG175$cd40)/ACTG175$cd40)[!ex]
summary(y)

## 0=failure, 1=success
y <- 1*(y > -0.5)

## summarize CD4 outcomes
table(y, ACTG175$arms[!ex])

table(y, ACTG175$arms[!ex])/
    matrix(table(ACTG175$arms[!ex]), nrow=2, ncol=4, byrow=TRUE)

## drop unneeded and unwanted variables
## 1: 'pidnum' patient ID number
##14: 'str2' which will be handled by strat1 below
##15: 'strat' which will be handled by strat1-strat3 below
##17: 'treat' handled by arm0-arm3 below
##18: 'offtrt' indicator of off-treatment before 96 weeks
##20: 'cd420' CD4 T cell count at 20 weeks
##21: 'cd496' CD4 T cell count at 96 weeks
##22: 'r' missing CD4 T cell count at 96 weeks
##24: 'cd820' CD8 T cell count at 20 weeks
##25: 'cens' indicator of observing the event in days
##26: 'days' number of days until the primary endpoint
##27: 'arms' handled by arm0-arm3 below
train <- as.matrix(ACTG175)[!ex, -c(1, 14:15, 17, 18, 20:22, 24:27)]
train <- cbind(1*(ACTG175$strat[!ex]==1), 1*(ACTG175$strat[!ex]==2),
               1*(ACTG175$strat[!ex]==3), train)
dimnames(train)[[2]][1:3] <- paste0('strat', 1:3)
train <- cbind(1*(ACTG175$arms[!ex]==0), 1*(ACTG175$arms[!ex]==1),
               1*(ACTG175$arms[!ex]==2), 1*(ACTG175$arms[!ex]==3), train)
dimnames(train)[[2]][1:4] <- paste0('arm', 0:3)

N <- nrow(train)

test0 <- train; test0[ , 1:4] <- 0; test0[ , 1] <- 1
test1 <- train; test1[ , 1:4] <- 0; test1[ , 2] <- 1
test2 <- train; test2[ , 1:4] <- 0; test2[ , 3] <- 1
test3 <- train; test3[ , 1:4] <- 0; test3[ , 4] <- 1

test <- rbind(test0, test1, test2, test3)

##test BART with token run to ensure installation works
## set.seed(21)
## post <- lbart(train, y, test, nskip=5, ndpost=5)

\dontrun{
set.seed(21)
post <- lbart(train, y, test)

## turn z-scores into probabilities
post$prob.test <- plogis(post$yhat.test)

## average over the posterior samples
post$prob.test.mean <- apply(post$prob.test, 2, mean)

## place estimates for arms 0-3 next to each other for convenience
itr <- cbind(post$prob.test.mean[(1:N)], post$prob.test.mean[N+(1:N)],
             post$prob.test.mean[2*N+(1:N)], post$prob.test.mean[3*N+(1:N)])

## find the BART ITR for each patient
itr.pick <- integer(N)
for(i in 1:N) itr.pick[i] <- which(itr[i, ]==max(itr[i, ]))-1

## arms 0 and 3 (monotherapy) are never chosen
table(itr.pick)

## do arms 1 and 2 show treatment heterogeneity?
diff. <- apply(post$prob.test[ , 2*N+(1:N)]-post$prob.test[ , N+(1:N)], 2, mean)
plot(sort(diff.), type='h', main='ACTG175 trial: 50\% CD4 decline from baseline at 96 weeks',
     xlab='Arm 2 (1) Preferable to the Right (Left)', ylab='Prob.Diff.: Arms 2 - 1')

library(rpart)
library(rpart.plot)

## make data frame for nicer names in the plot
var <- as.data.frame(train[ , -(1:4)])

dss <- rpart(diff. ~ var$age+var$gender+var$race+var$wtkg+var$cd40+var$cd80+
                   var$karnof+var$symptom+var$hemo+var$homo+var$drugs+var$z30+
                   var$zprior+var$oprior+var$strat1+var$strat2+var$strat3,
               method='anova', control=rpart.control(cp=0.1))
rpart.plot(dss, type=3, extra=101)

## if strat1==1 (antiretroviral naive), then arm 2 is better
## otherwise, arm 1
print(dss)

all0 <- apply(post$prob.test[ , (1:N)], 1, mean)
all1 <- apply(post$prob.test[ , N+(1:N)], 1, mean)
all2 <- apply(post$prob.test[ , 2*N+(1:N)], 1, mean)
all3 <- apply(post$prob.test[ , 3*N+(1:N)], 1, mean)

## BART ITR
BART.itr <- apply(post$prob.test[ , c(N+which(itr.pick==1), 2*N+which(itr.pick==2))], 1, mean)

test <- train
test[ , 1:4] <- 0
test[test[ , 5]==0, 2] <- 1
test[test[ , 5]==1, 3] <- 1

## BART ITR simple
BART.itr.simp <- pwbart(test, post$treedraws)
BART.itr.simp <- apply(plogis(BART.itr.simp), 1, mean)

plot(density(BART.itr), xlab='Value', xlim=c(0.475, 0.775), lwd=2,
     main='ACTG175 trial: 50\% CD4 decline from baseline at 96 weeks')
lines(density(BART.itr.simp), col='brown', lwd=2)
lines(density(all0), col='green', lwd=2)
lines(density(all1), col='red', lwd=2)
lines(density(all2), col='blue', lwd=2)
lines(density(all3), col='yellow', lwd=2)
legend('topleft', legend=c('All Arm 0 (ZDV only)', 'All Arm 1 (ZDV+DDI)',
                           'All Arm 2 (ZDV+DDC)', 'All Arm 3 (DDI only)',
                           'BART ITR simple', 'BART ITR'),
       col=c('green', 'red', 'blue', 'yellow', 'brown', 'black'), lty=1, lwd=2)

}
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
