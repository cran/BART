\name{mc.wbart.gse}
\title{Global SE variable selection for BART with parallel computation}
\alias{mc.wbart.gse}
\description{
Here we implement the global SE method for variable selection in
nonparametric survival analysis with BART.  Unfortunately, the method is
very computationally intensive so we present some trade-offs below. }

\usage{
mc.wbart.gse( x.train, y.train,
              P=50L, R=5L, ntree=20L, numcut=100L, C=1, alpha=0.05,
              k=2.0, power=2.0, base=0.95,
              ndpost=2000L, nskip=100L,
              printevery=100L, keepevery=1L, keeptrainfits=FALSE,
              seed=99L, mc.cores=2L, nice=19L 
              )
}

\arguments{

  \item{x.train}{ Explanatory variables for training (in sample)
data.\cr Must be a matrix with (as usual) rows corresponding to
observations and columns to variables.\cr \code{surv.bart} will generate
draws of \eqn{f(t, x)} for each \eqn{x} which is a row of x.train.  }

   \item{y.train}{
     The continuous outcome.
   }

   \item{P}{
     The number of permutations: typically 50 or 100.
   }

   \item{R}{
     The number of replicates: typically 5 or 10.
   }   

   \item{ntree}{
     The number of trees.  In variable selection,
     the number of trees is smaller than what might
     be used for the best fit.
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

   \item{C}{
     The starting value for the multiple of SE.  You should not need to
     change this except in rare circumstances.
   }
   
   \item{alpha}{
     The global SE method relies on simultaneous 1-\code{alpha} coverage
     across the permutations for all predictor variables.
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

   \item{ndpost}{
   The number of posterior draws after burn in.  In the global SE
   method, generally, the method is repeated several times to
   establish the variable count probabilities.  However, we take the
   alternative approach of simply running the MCMC chain longer which
   should result in the same stabilization of the estimates.  Therefore,
   the number of posterior draws in variable selection should be set to a larger value than
   would be typically anticipated for fitting. 
   }

   \item{nskip}{
   Number of MCMC iterations to be treated as burn in.
   }

   \item{printevery}{
   As the MCMC runs, a message is printed every printevery draws.
   }

   \item{keepevery}{
   Every \code{keepevery} draw is kept.
   }

   \item{keeptrainfits}{
   If \code{TRUE} the draws of \eqn{f(t, x)} for \eqn{x} = rows of
   x.train are generated.
   }

   \item{seed}{
     \code{seed} required for reproducible MCMC.
   }

   \item{mc.cores}{
     Number of cores to employ in parallel.
   }
   
   \item{nice}{
     Set the job priority.  The default
     priority is 19: priorities go from 0 (highest) to 19 (lowest).
   }
}

\value{
   \code{mc.wbart.gse} returns a list.
}
\references{
  Bleich, J., Kapelner, A., George, E.I., and Jensen, S.T. (2014).
  Variable selection for BART: an application to gene regulation.
  \emph{The Annals of Applied Statistics}, \bold{8:1750-81}.
}
%% \author{
%% Rodney Sparapani: \email{rsparapa@mcw.edu}
%% }
\seealso{
\code{\link{mc.wbart}}
}
\examples{
\dontrun{

library(ElemStatLearn)

data(phoneme)

x.train <- matrix(NA, nrow=4509, ncol=257)
    
dimnames(x.train)[[2]] <- c(paste0('x.', 1:256), 'speaker')
    
x.train[ , 257] <- as.numeric(phoneme$speaker)

for(j in 1:256) x.train[ , j] <- as.numeric(phoneme[ , paste0('x.', j)])

gse <- mc.wbart.gse(x.train, as.numeric(phoneme$g), mc.cores=5, seed=99)

## important variables
dimnames(x.train)[[2]][gse$which]
    
}
}
\keyword{nonproportional hazards variable selection}
