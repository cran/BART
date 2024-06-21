
% copyright (C) 1994-1999 Bill Venables and Brian Ripley
%
\name{class.ind}
\alias{class.ind}
\title{
Generates Class Indicator Matrix from a Factor
}
\description{
Generates a class indicator function from a given factor.
}
\usage{
class.ind(cl)
}
\arguments{
\item{cl}{
factor or vector of classes for cases.
}}
\value{
a matrix which is zero except for the column corresponding to the class.
}
\references{
  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S.} Fourth edition.  Springer.
}
\examples{

set.seed(99)

a <- rbinom(20, 4, 0.5)

table(a)

b <- class.ind(a)

str(b)

t(cbind(a, b))

}
\keyword{utilities}
