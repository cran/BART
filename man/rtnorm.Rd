\name{rtnorm}
\title{Testing truncated Normal sampling}
\alias{rtnorm}
%\alias{rtnorm_reject}

\description{ Truncated Normal latents are necessary to transform a
  binary BART into a continuous BART.  
}

\usage{
rtnorm(mean, tau, sd)
%rtnorm_reject(mean, tau, sd)
}

\arguments{
   \item{mean}{ Mean of truncated Normal. }
   \item{tau}{ Truncation point of truncated Normal. }
   \item{sd}{ Standard deviation of truncated Normal. }
 }
 
\value{
  Returns a truncated Normal, i.e., \eqn{N(mean, sd)I(tau, infinity)}.
}

\references{
Robert, C.P. (1995)
   Simulation of truncated normal variables.
   \emph{Statistics and computing}, \bold{5(2)}, 121--125.
}

\author{
Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
Rodney Sparapani: \email{rsparapa@mcw.edu}\cr,
Robert Gramacy: \email{rbg@vt.edu}.
}
\seealso{
\code{\link{pbart}, \link{lbart}}
}
\examples{

set.seed(12)

rtnorm(0, 3, 1)
rtnorm(0, 3, 1)


}

