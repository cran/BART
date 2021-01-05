\name{draw_lambda_i}
\title{Testing truncated Normal sampling}
\alias{draw_lambda_i}
%\alias{draw_lambda_i_reject}

\description{ Truncated Normal latents with non-unit variance
  are necessary for logistic BART.  }

\usage{
draw_lambda_i(lambda, mean, kmax=1000, thin=1)
%draw_lambda_i_reject(mean, tau, sd)
}

\arguments{
   \item{lambda}{ Previous value of lambda.}
   \item{mean}{ Mean of truncated Normal. }
   \item{kmax}{ The number of terms in the mixture. }
   \item{thin}{ The thinning parameter. }
 }
 
\value{ Returns the variance for a truncated Normal, i.e., \eqn{N(mean,
  lambda)I(tau, infinity)}.  }

%% \author{
%% Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
%% Rodney Sparapani: \email{rsparapa@mcw.edu}\cr,
%% Robert Gramacy: \email{rbg@vt.edu}.
%% }
\seealso{
\code{\link{rtnorm}, \link{lbart}}
}
\examples{

set.seed(12)

draw_lambda_i(1, 2)
rtnorm(1, 2, sqrt(6.773462), 6)
draw_lambda_i(6.773462, 2)

}

