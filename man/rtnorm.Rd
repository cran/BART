\name{rtnorm}
\title{Testing truncated Normal sampling}
\alias{rtnorm}
%\alias{rtnorm_reject}

\description{ Truncated Normal latents are necessary to transform a
  binary BART into a continuous BART.  
}

\usage{
rtnorm(n, mean, sd, tau)
%rtnorm_reject(mean, tau, sd)
}

\arguments{
   \item{n}{ Number of samples. }
   \item{mean}{ Mean. }
   \item{sd}{ Standard deviation. }
   \item{tau}{ Truncation point. }
 }
 
\value{
  Returns \code{n} truncated Normals, i.e., \eqn{N(mean, sd)I(tau, infinity)}.
}

\references{
  Robert C. (1995)
  Simulation of truncated normal variables.
  \emph{Statistics and computing}, \bold{5(2)}, 121--125.
}

%% \author{
%% Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
%% Rodney Sparapani: \email{rsparapa@mcw.edu}\cr,
%% Robert Gramacy: \email{rbg@vt.edu}.
%% }
\seealso{
\code{\link{pbart}, \link{lbart}}
}
\examples{

set.seed(12)

rtnorm(1, 0, 1, 3)
rtnorm(1, 0, 1, 3)


}

