\name{rtgamma}
\title{Testing truncated Gamma sampling}
\alias{rtgamma}

\description{ Truncated Gamma draws are needed for the standard
  deviation of the random effects Gibbs conditional.
}

\usage{
rtgamma(shape, rate, a)
}

\arguments{
  \item{shape}{ Sampling from a truncated Gamma where
    \eqn{E[x]=shape/rate}.}
  \item{rate}{ This parameter is the inverse of the scale
    which is an alternative representation for the Gamma distribution. }
   \item{a}{ The truncation point, i.e., \eqn{a<x}. }
 }
 
\value{
  Returns a truncated Gamma, i.e., \eqn{Gam(shape, rate)I(a, infinity)}.
}

\references{
Gentle, J. E. (2013)
   Random number generation and Monte Carlo methods. Springer, New York, NY.
}

\author{
Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
Rodney Sparapani: \email{rsparapa@mcw.edu}\cr,
Robert Gramacy: \email{rbg@vt.edu}.
}

%% \seealso{
%% \code{\link{gbmm}}
%% }
\examples{

set.seed(12)
rtgamma(3, 1, 4)
rtgamma(3, 1, 4)

N=10000
a=double(N)
for(i in 1:N) a[i]=rtgamma(10, 2, 1)
mean(a)
min(a)

}

