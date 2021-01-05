\name{crisk2.bart}
\title{BART for competing risks}
\alias{crisk2.bart}
\alias{mc.crisk2.bart}
\description{
Here we have implemented another approach to utilize BART
for competing risks that is very flexible, and is akin to discrete-time
survival analysis.  Following the capabilities of BART, we allow for
maximum flexibility in modeling the dependence of competing failure times on
covariates.  In particular, we do not impose proportional hazards.

Similar to \code{crisk.bart}, we utilize two BART models, yet they are
two different BART models than previously considered.  First, given an
event of either cause occurred, we employ a typical binary BART model to
discriminate between cause 1 and 2.  Next, we proceed as if it were a typical
survival analysis with BART for an absorbing event from either cause.

To elaborate, consider data in the form: \eqn{(s_i, \delta_i,
{x}_i)}{(s, delta, x)} where \eqn{s_i}{s} is the event time;
\eqn{\delta_i}{delta} is an indicator distinguishing events,
\eqn{\delta_i=h}{delta=h} due to cause \eqn{h in {1, 2}}, from
right-censoring, \eqn{\delta_i=0}{delta=0}; \eqn{{x}_i}{x} is a vector of
covariates; and \eqn{i=1, ..., N}{i=1, ..., N (i suppressed for
convenience)} indexes subjects.
We denote the \eqn{K} distinct event/censoring times by
\eqn{0<t_{(1)}<...<t_{(K)}<\infty}{0<t(1)<...< t(K)<infinity} thus
taking \eqn{t_{(j)}}{t(j)} to be the \eqn{j^{th}}{j'th} order statistic
among distinct observation times and, for convenience,
\eqn{t_{(0)}=0}{t(0)=0}. 
%\eqn{n_i=\sum_j I[t_{(j)}\leq t_i]}{n=sum I[t(j)<=t]}
%% This means \eqn{y_{ij}=0}{y(j)=0} if
%% \eqn{j<n_i}{j<n} and \eqn{y_{in_i}=\delta_i}{y(n)=delta}.

First, consider event indicators for an event from either cause:
\eqn{y_{1ij}}{y(1j)} for each subject \eqn{i} at each distinct time
\eqn{t_{(j)}}{t(j)} up to and including the subject's last observation
time \eqn{s_i=t_{(n_i)}}{s=t(n)} with \eqn{n_i=\arg \max_j [t_{(j)}\leq
s_i]}{n=arg max [t(j)<=s]}.  We denote by \eqn{p_{1ij}}{p(1j)} the
probability of an event at time \eqn{t_{(j)}}{t(j)} conditional on no
previous event. We now write the model for \eqn{y_{1ij}}{y(1j)} as a
nonparametric probit (or logistic) regression of \eqn{y_{1ij}}{y(1j)} on
the time \eqn{t_{(j)}}{t(j)} and the covariates \eqn{{x}_{1i}}{x(1)},
and then utilize BART for binary responses.  Specifically, \eqn{
y_{1ij}\ =\ I[\delta_i>0] I[s_i=t_{(j)}],\ j=1, ..., n_i}{ y(1j) =
I[delta>0] I[s=t(j)], j=1, ..., n}.  Therefore, we have \eqn{p_{1ij} =
F(mu_{1ij}),\ mu_{1ij} = mu_1+f_1(t_{(j)}, {x}_{1i})}{p(1j) = F(mu(1j)),
mu(1j) = mu(h)+f1(t(j), x(1))} where \eqn{F} denotes the Normal (or
Logistic) cdf.  

Next, we denote by \eqn{p_{2i}}{p(2)} the probability of a cause 1
event at time \eqn{s_i}{s} conditional on an event having
occurred.  We now write the model for \eqn{y_{2i}}{y(2)} as a
nonparametric probit (or logistic) regression of \eqn{y_{2i}}{y(2)} on
the time \eqn{s_i}{s} and the covariates \eqn{{x}_{2i}}{x(2)},
via BART for binary responses.  Specifically, \eqn{
y_{2i}\ =\ I[\delta_i=1]}{ y(2) = I[delta=1]}.  Therefore, we
have \eqn{p_{2i} = F(mu_{2i}),\ mu_{2i} = mu_2+f_2(s_i,
{x}_{2i})}{p(2j) = F(mu(2j)), mu(2j) = mu(2)+f2(s, x(h))} where
\eqn{F} denotes the Normal (or Logistic) cdf.  Although, we modeled
\eqn{p_{2i}}{p(2)} at the time of an event, \eqn{s_i}{s}, we can
estimate this probability at any other time points on the grid via
\eqn{p(t_{(j)}, x_2)=F( mu_2+f_2(t_{(j)}, {x}_2))}{p(t(j), x(2)) = F(
  mu(2)+f2(t(j), x(2)))}.
Finally, based on these probabilities,
\eqn{p_{hij}}{p(hj)}, we can construct targets of inference such as the
cumulative incidence functions.
 }

\usage{

crisk2.bart(x.train=matrix(0,0,0), y.train=NULL,
           x.train2=x.train, y.train2=NULL,
           times=NULL, delta=NULL, K=NULL,
           x.test=matrix(0,0,0), x.test2=x.test, 
           sparse=FALSE, theta=0, omega=1,
           a=0.5, b=1, augment=FALSE,
           rho=NULL, rho2=NULL,
           xinfo=matrix(0,0,0), xinfo2=matrix(0,0,0),
           usequants=FALSE, %cont=FALSE,
           rm.const=TRUE, type='pbart',
           ntype=as.integer(
               factor(type, levels=c('wbart', 'pbart', 'lbart'))),
           k=2, power=2, base=0.95,
           offset=NULL, offset2=NULL,
           tau.num=c(NA, 3, 6)[ntype],
           %binaryOffset=NULL, binaryOffset2=NULL,
           ntree=50, numcut=100, ndpost=1000, nskip=250,
           keepevery = 10L,
           %nkeeptrain=ndpost,
           %nkeeptest=ndpost,
           %nkeeptestmean=ndpost,
           %nkeeptreedraws=ndpost,
           printevery=100L, %treesaslists=FALSE,
           %keeptrainfits=TRUE,
           id=NULL,    ## crisk2.bart only
           seed=99,    ## mc.crisk2.bart only
           mc.cores=2, ## mc.crisk2.bart only
           nice=19L    ## mc.crisk2.bart only
          )

mc.crisk2.bart(x.train=matrix(0,0,0), y.train=NULL,
              x.train2=x.train, y.train2=NULL,
              times=NULL, delta=NULL, K=NULL,
              x.test=matrix(0,0,0), x.test2=x.test, 
              sparse=FALSE, theta=0, omega=1,
              a=0.5, b=1, augment=FALSE,
              rho=NULL, rho2=NULL,
              xinfo=matrix(0,0,0), xinfo2=matrix(0,0,0),
              usequants=FALSE, %cont=FALSE,
              rm.const=TRUE, type='pbart',
              ntype=as.integer(
                  factor(type, levels=c('wbart', 'pbart', 'lbart'))),
              k=2, power=2, base=0.95,
              offset=NULL, offset2=NULL,
              tau.num=c(NA, 3, 6)[ntype],
              %binaryOffset=NULL, binaryOffset2=NULL,
              ntree=50, numcut=100, ndpost=1000, nskip=250,
              keepevery = 10L,
              %nkeeptrain=ndpost, nkeeptest=ndpost,
              %nkeeptestmean=ndpost,
              %nkeeptreedraws=ndpost,
              printevery=100L, %treesaslists=FALSE,
              %keeptrainfits=TRUE,
              id=NULL,    ## crisk2.bart only
              seed=99,    ## mc.crisk2.bart only
              mc.cores=2, ## mc.crisk2.bart only
              nice=19L    ## mc.crisk2.bart only
             )
}

\arguments{

  \item{x.train}{ Covariates for training (in sample) data for an event.\cr
    Must be a data.frame or a matrix with rows corresponding to
observations and columns to variables.\cr \code{crisk2.bart} will
generate draws of \eqn{f_1(t, x)}{f1(t, x)} for each \eqn{x} which is
a row of \code{x.train} (note that the definition of \code{x.train} is
dependent on whether \code{y.train} has been specified; see below).  }

   \item{y.train}{ Event binary response for training (in sample)
   data.\cr If \code{y.train} is \code{NULL}, then \code{y.train}
   (\code{x.train} and \code{x.test}, if specified) are generated by a
   call to \code{surv.pre.bart} (which require that \code{times} and
   \code{delta} be provided: see below); otherwise, \code{y.train}
   (\code{x.train} and \code{x.test}, if specified) are utilized as
   given assuming that the data construction has already been performed.
   }

   \item{x.train2}{ Covariates for training (in sample)
     data of for a cause 1 event. Similar to \code{x.train} above.}

   \item{y.train2}{
   Cause 1 event binary response for training (in sample) data.
   Similar to \code{y.train} above.
   }

   \item{times}{
     The time of event or right-censoring, \eqn{s_i}{s}.\cr
     If \code{y.train} is \code{NULL}, then \code{times} (and \code{delta})
     must be provided.
   }

   \item{delta}{
     The event indicator: 1 for cause 1, 2 for cause 2 and 0 is censored.\cr
     If \code{y.train} is \code{NULL}, then \code{delta} (and \code{times})
     must be provided.
   }

   \item{K}{
     If provided, then coarsen \code{times} per the quantiles
     \eqn{1/K, 2/K, ..., K/K}.
   }
   
   \item{x.test}{
   Covariates for test (out of sample) data of an event.\cr
   Must be a data.frame or a matrix and have the same structure as
   \code{x.train}.\cr
   \code{crisk2.bart} will generate draws of
   \eqn{f_1(t, x)}{f1(t, x)} for each \eqn{x} which is a row of
   \code{x.test}.
 }

    \item{x.test2}{
      Covariates for test (out of sample) data of a cause 1 event.
     Similar to \code{x.test} above.}

   \item{sparse}{Whether to perform variable selection based on a
     sparse Dirichlet prior; see Linero 2016.}
   \item{theta}{Set \eqn{theta} parameter; zero means random.}
   \item{omega}{Set \eqn{omega} parameter; zero means random.}
   \item{a}{Sparse parameter for \eqn{Beta(a, b)} prior:
     \eqn{0.5<=a<=1} where lower values inducing more sparsity.}
   \item{b}{Sparse parameter for \eqn{Beta(a, b)} prior; typically,
     \code{b=1}.}
   \item{rho}{Sparse parameter: typically \code{rho=p} where \code{p} is the
     number of covariates in \code{x.train}.}
   \item{rho2}{Sparse parameter: typically \code{rho2=p} where \code{p} is the
     number of covariates  in \code{x.train2}.}
   \item{augment}{Whether data augmentation is to be performed in sparse
     variable selection.}
   
   \item{xinfo}{ You can provide the cutpoints to BART or let BART
     choose them for you.  To provide them, use the \code{xinfo}
     argument to specify a list (matrix) where the items (rows) are the
     covariates and the contents of the items (columns) are the
     cutpoints.  }

   \item{xinfo2}{ Cause 2 cutpoints. }

   \item{usequants}{ If \code{usequants=FALSE}, then the
    cutpoints in \code{xinfo} are generated uniformly; otherwise,
    if \code{TRUE}, uniform quantiles are used for the cutpoints. }
   
   %\item{cont}{ Whether or not to assume all variables are continuous.}

   \item{rm.const}{ Whether or not to remove constant variables.}

   \item{type}{ Whether to employ probit BART via Albert-Chib,
     \code{'pbart'}, or logistic BART by Holmes-Held, \code{'lbart'}. }
 \item{ntype}{ The integer equivalent of \code{type} where
  \code{'wbart'} is 1, \code{'pbart'} is 2 and
  \code{'lbart'} is 3.}
         
   \item{k}{ k is the number of prior standard deviations
     \eqn{f_h(t, x)}{f(h)(t, x)}
   is away from +/-3.  The bigger k is, the more conservative the
   fitting will be.  }

   \item{power}{
   Power parameter for tree prior.
   }

   \item{base}{
   Base parameter for tree prior.
   }

   \item{offset}{
   Cause 1 binary offset.
   }

   \item{offset2}{
   Cause 2 binary offset.
   }
   \item{tau.num}{ The numerator in the \code{tau} definition. }
   \item{ntree}{
   The number of trees in the sum.
   }

   \item{numcut}{ The number of possible values of cutpoints (see
   \code{usequants}).  If a single number if given, this is used for all
   variables.  Otherwise a vector with length equal to
   \code{ncol(x.train)} is required, where the \eqn{i^{th}}{i'th}
   element gives the number of cutpoints used for the \eqn{i^{th}}{i'th}
   variable in \code{x.train}.  If \code{usequants} is \code{FALSE},
   \code{numcut} equally spaced cutoffs are used covering the range of
   values in the corresponding column of \code{x.train}.  If
   \code{usequants} is \code{TRUE}, then \code{min(numcut, the number of
   unique values in the corresponding columns of x.train - 1)} cutpoint
   values are used.  }

   \item{ndpost}{
   The number of posterior draws returned.
   }

   \item{nskip}{
   Number of MCMC iterations to be treated as burn in.
   }

   \item{keepevery}{
   Every \code{keepevery} draw is kept to be returned to the user.
   }

   %% \item{nkeeptrain}{
   %% Number of MCMC iterations to be returned for training data.
   %% }

   %% \item{nkeeptest}{
   %% Number of MCMC iterations to be returned for test data.
   %% }

   %% \item{nkeeptestmean}{
   %% Number of MCMC iterations to be returned for test mean.
   %% }

   %% \item{nkeeptreedraws}{
   %% Number of MCMC iterations to be returned for tree draws.
   %% }

   \item{printevery}{
   As the MCMC runs, a message is printed every \code{printevery} draws.
   }

   %% \item{keeptrainfits}{
   %% Whether to keep \code{yhat.train} or not.
   %% }

   %% \item{treesaslists}{
   %%   The default behavior is to return the trees as a text string.  If you
   %%   specify \code{treesaslists=TRUE}, then the trees will be returned
   %%   as a list as well.
   %% }

   \item{id}{
     \code{crisk2.bart} only: unique identifier added to returned list.
   }

   \item{seed}{
     \code{mc.crisk2.bart} only: seed required for reproducible MCMC.
   }

   \item{mc.cores}{
     \code{mc.crisk2.bart} only: number of cores to employ in parallel.
   }

   \item{nice}{ \code{mc.crisk2.bart} only: set the job niceness.  The
     default niceness is 19: niceness goes from 0 (highest priority) to
     19 (lowest priority).  }
 }

\value{
  \code{crisk2.bart} returns an object of type \code{crisk2bart} which is
  essentially a list.  Besides the items listed
   below, the list has \code{offset}, \code{offset2},
   \code{times} which are the unique times, \code{K}
   which is the number of unique times, \code{tx.train} and
  \code{tx.test}, if any.

   \item{yhat.train}{ A matrix with \code{ndpost} rows and
   \code{nrow(x.train)} columns.  Each row corresponds to a draw
   \eqn{f^*_1}{f1*} from the posterior of \eqn{f_1}{f1} and each column
   corresponds to a row of \code{x.train}.  The \eqn{(i,j)} value is
   \eqn{f^*_1(t, x)}{f1*(t, x)} for the \eqn{i^{th}}{i'th} kept draw of
   \eqn{f_1}{f1} and the \eqn{j^{th}}{j'th} row of \code{x.train}.
   Burn-in is dropped.  }

   \item{yhat.test}{Same as \code{yhat.train} but now the x's are the rows of the test data.}
   \item{surv.test}{test data fits for the survival function, \eqn{S(t, x)}.}
   \item{surv.test.mean}{mean of \code{surv.test} over the posterior samples.}
   \item{prob.test}{The probability of suffering an event.}
   \item{prob.test2}{The probability of suffering a cause 1 event.}
   \item{cif.test}{The cumulative incidence function of cause 1,
   \eqn{F_1(t, x)}{F1(t, x)}.}
   \item{cif.test2}{The cumulative incidence function of cause 2,
   \eqn{F_2(t, x)}{F2(t, x)}.}
   %\item{yhat.train.mean}{train data fits = mean of yhat.train columns.}
   %\item{yhat.test.mean}{test data fits = mean of yhat.test columns.}
   \item{cif.test.mean}{mean of \code{cif.test} columns for cause 1.}
   \item{cif.test2.mean}{mean of \code{cif.test2} columns for cause 2.}
   \item{varcount}{a matrix with \code{ndpost} rows and
   \code{nrow(x.train)} columns.
   Each row is for a draw. For each variable (corresponding to the columns),
   the total count of the number of times
   this variable is used for an event in a tree decision rule (over all trees) is given.}
   \item{varcount2}{ For each variable the total count of the number of times
   this variable is used for a cause 1 event in a tree decision rule is given.}
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

%% Sparapani, Logan, McCulloch, Laud (2018)
%%   Nonparametric competing risks analysis using
%%   Bayesian Additive Regression Trees (BART).
%%   arXiv preprint arXiv:1806.11237

%% }
%% \author{
%% Rodney Sparapani: \email{rsparapa@mcw.edu}
%% }
\seealso{
\code{\link{surv.pre.bart}}, \code{\link{predict.crisk2bart}},
   \code{\link{mc.crisk2.pwbart}}, \code{\link{crisk.bart}}
}
\examples{

data(transplant)

pfit <- survfit(Surv(futime, event) ~ abo, transplant)

# competing risks for type O
plot(pfit[4,], xscale=7, xmax=735, col=1:3, lwd=2, ylim=c(0, 1),
       xlab='t (weeks)', ylab='Aalen-Johansen (AJ) CI(t)')
    legend(450, .4, c("Death", "Transplant", "Withdrawal"), col=1:3, lwd=2)
## plot(pfit[4,], xscale=30.5, xmax=735, col=1:3, lwd=2, ylim=c(0, 1),
##        xlab='t (months)', ylab='Aalen-Johansen (AJ) CI(t)')
##     legend(450, .4, c("Death", "Transplant", "Withdrawal"), col=1:3, lwd=2)

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

##test BART with token run to ensure installation works
set.seed(99)
post <- crisk2.bart(x.train=x.train, times=times, delta=delta,
                   x.test=x.test, nskip=1, ndpost=1, keepevery=1)

\dontrun{

## run one long MCMC chain in one process
## set.seed(99)
## post <- crisk2.bart(x.train=x.train, times=times, delta=delta, x.test=x.test)

## in the interest of time, consider speeding it up by parallel processing
## run "mc.cores" number of shorter MCMC chains in parallel processes
post <- mc.crisk2.bart(x.train=x.train, times=times, delta=delta,
                      x.test=x.test, seed=99, mc.cores=8)

K <- post$K

typeO.cif.mean <- apply(post$cif.test, 2, mean)
typeO.cif.025 <- apply(post$cif.test, 2, quantile, probs=0.025)
typeO.cif.975 <- apply(post$cif.test, 2, quantile, probs=0.975)

plot(pfit[4,], xscale=7, xmax=735, col=1:3, lwd=2, ylim=c(0, 0.8),
       xlab='t (weeks)', ylab='CI(t)')
points(c(0, post$times)*7, c(0, typeO.cif.mean), col=4, type='s', lwd=2)
points(c(0, post$times)*7, c(0, typeO.cif.025), col=4, type='s', lwd=2, lty=2)
points(c(0, post$times)*7, c(0, typeO.cif.975), col=4, type='s', lwd=2, lty=2)
     legend(450, .4, c("Transplant(BART)", "Transplant(AJ)",
                       "Death(AJ)", "Withdrawal(AJ)"),
            col=c(4, 2, 1, 3), lwd=2)
##dev.copy2pdf(file='../vignettes/figures/liver-BART.pdf')
## plot(pfit[4,], xscale=30.5, xmax=735, col=1:3, lwd=2, ylim=c(0, 0.8),
##        xlab='t (months)', ylab='CI(t)')
## points(c(0, post$times)*30.5, c(0, typeO.cif.mean), col=4, type='s', lwd=2)
## points(c(0, post$times)*30.5, c(0, typeO.cif.025), col=4, type='s', lwd=2, lty=2)
## points(c(0, post$times)*30.5, c(0, typeO.cif.975), col=4, type='s', lwd=2, lty=2)
##      legend(450, .4, c("Transplant(BART)", "Transplant(AJ)",
##                        "Death(AJ)", "Withdrawal(AJ)"),
##             col=c(4, 2, 1, 3), lwd=2)

}
}
\keyword{nonparametric survival model}
\keyword{nonproportional hazards}
