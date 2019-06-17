\name{stratrs}
\alias{stratrs}
\title{Perform stratified random sampling to balance outcomes}
\description{This function is used to perform stratified random
  sampling to balance outcomes among the shards.
}
\usage{
stratrs(y, C=5, P=0)
}

\arguments{
  \item{y}{The binary/categorical/continuous outcome.}
  \item{C}{The number of shards to break the data set into.}
  \item{P}{For continuous data, we break the range into
    P segments via the quantiles.  Specifying, P=20 seems to
    work reasonably well.}
}

\details{
To perform BART with large data sets, random sampling is employed
  to break the data into \code{C} shards.  Each shard should be
  balanced with respect to the outcome.  For binary/categorical
  outcomes, stratified random sampling is employed with this function.
}
\value{
A vector is returned with each element assigned to a shard.
}

\seealso{
   \code{\link{rs.pbart}}
 }

 \examples{
set.seed(12)
x <- rbinom(25000, 1, 0.1)
a <- stratrs(x)
table(a, x)
z <- pmin(rpois(25000, 0.8), 5)
b <- stratrs(z)
table(b, z)
}
\keyword{convergence diagnostics}
