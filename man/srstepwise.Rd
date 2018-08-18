\name{srstepwise}
\alias{srstepwise}
\title{Stepwise Variable Selection Procedure for survreg}
\usage{
srstepwise(x, times, delta, sle = 0.15, sls = 0.15, dist='lognormal')
}
\arguments{
\item{x}{Matrix of variables to consider.}
  
\item{times}{The time to an event, if any.}

\item{delta}{The event indicator: 1 for event, 0 for no event.}

\item{sle}{The chosen significance level for entering.}

\item{sls}{The chosen significance level for staying.}

\item{dist}{The distribution to be used by \code{survreg}.}
}
\value{
Returns a list of indices of variables which have entered and stayed.
}
\description{
This stepwise variable selection procedure can be applied to obtain the
best candidates for a \code{survreg} fit.
}
\details{
  Unfortunately, no stepwise procedure exists for \code{survreg} models.
  Therefore, we provide this brute force method.
}

\examples{

names. <- names(lung)[-(2:3)]
status1 <- ifelse(lung$status==2,1,0)
X <- as.matrix(lung)[ , names.]
vars=srstepwise(X, lung$time, status1)
print(names.[vars])

}
\seealso{
\link{lung}
}
