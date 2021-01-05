\name{mc.cores.openmp}
\title{Detecting OpenMP}
\alias{mc.cores.openmp}

\description{ This package was designed for OpenMP.  For example, the
  \code{pwbart} function can use OpenMP or the parallel R package for
  multi-threading.  On UNIX/Unix-like systems, OpenMP, if available, is
  discovered at install time; for the details, see the
  \code{configure.ac} file which can be found in the source version of
  this package.  However, we know of no GPL licensed code available to
  detect OpenMP on Windows (for Artistic licensed OpenMP detection code
  on Windows, see the Bioconductor R package rGADEM).  To determine
  whether OpenMP is available at run time, we provide the function
  documented here.  }

\usage{
mc.cores.openmp()
}

\value{
  Returns a zero when OpenMP is not available, otherwise, an integer
  greater than zero when OpenMP is available (returns one unless
  you are running in a multi-threaded process).
}

%% \author{
%% Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
%% Rodney Sparapani: \email{rsparapa@mcw.edu}.
%% }
\seealso{
\code{\link{pwbart}}
}
\examples{

mc.cores.openmp()

}

\keyword{multi-threading}
\keyword{OpenMP}
\keyword{parallel}
