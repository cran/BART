\name{spectrum0ar}
\alias{spectrum0ar}
\title{Estimate spectral density at zero}
\description{
The spectral density at frequency zero is estimated by fitting an
autoregressive model.  \code{spectrum0(x)/length(x)} estimates the
variance of \code{mean(x)}.
}
\usage{
spectrum0ar(x)
}

\arguments{ \item{x}{Matrix of MCMC chains: the rows are the samples and
  the columns are different "parameters".  For BART, generally, the
  columns are estimates of \eqn{f}.  For \code{pbart}, they are
  different subjects.  For \code{surv.bart}, they are different subjects
  at a grid of times.
}}

\details{
The \code{ar()} function to fit an autoregressive model to the time
series x. For multivariate time series, separate models are fitted for
each column. The value of the spectral density at zero is then given
by a well-known formula.  Adapted from the \code{spectrum0.ar} function of
the coda package which passes \code{mcmc} objects as arguments
rather than matrices.
}
\value{
A list with the following values
\item{spec}{The predicted value of the spectral density at frequency zero.}
\item{order}{The order of the fitted model}
}

\references{
Martyn Plummer, Nicky Best, Kate Cowles and Karen Vines (2006). 
CODA: Convergence Diagnosis and Output Analysis for MCMC, R News, vol 6,
7-11.

BW Silverman (1986).
Density estimation for statistics and data analysis.
Chapman and Hall, London.
}

\seealso{
   \code{\link{gewekediag}}
}
\keyword{convergence diagnostics}
