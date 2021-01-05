\name{rtgamma}
\title{Testing truncated Gamma sampling}
\alias{rtgamma}

\description{ Truncated Gamma draws are needed for the standard
  deviation of the random effects Gibbs conditional.
}

\usage{
rtgamma(n, shape, rate, a)
}

\arguments{
   \item{n}{ Number of samples. }
  \item{shape}{ Sampling from a truncated Gamma where
    \eqn{E[x]=shape/rate}.}
  \item{rate}{ This parameter is the inverse of the scale
    which is an alternative representation for the Gamma distribution. }
   \item{a}{ The truncation point, i.e., \eqn{a<x}. }
 }
 
\value{
  Returns \code{n} truncated Gamma, i.e., \eqn{Gam(shape, rate)I(a, infinity)}.
}

\references{
  Gentle J. (2013)
  Random number generation and Monte Carlo methods.
  Springer, New York, NY.
}

%% \author{
%% Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
%% Rodney Sparapani: \email{rsparapa@mcw.edu}\cr,
%% Robert Gramacy: \email{rbg@vt.edu}.
%% }

%% \seealso{
%% \code{\link{gbmm}}
%% }
\examples{

set.seed(12)
rtgamma(1, 3, 1, 4)
rtgamma(1, 3, 1, 4)

a=rtgamma(10000, 10, 2, 1)
mean(a)
min(a)

}

