\name{BART-package}
\alias{BART-package}
\alias{BART}
\docType{package}
\title{
  Bayesian Additive Regression Trees (BART) provide flexible nonparametric 
  modeling of covariates for continuous, binary, categorical and time-to-event outcomes.
}
\description{
  Efficient and reliable software is critical to the successful usage of the BART
  methodology.  We believe that this package fills that need.  The BART
  engine provided here is written in C++ for reasons of speed and
  maintainability; the source code is easy to read and modify while
  still being relatively efficient computationally.  This package brings
  together R functions for a variety of outcomes: continuous,
  dichotomous, categorical and time-to-event with right censoring.
}
\details{
For continuous outcomes, use the \code{wbart}/\code{mc.wbart} function(s). 
For dichotomous outcomes, use the \code{pbart}/\code{mc.pbart} or
\code{lbart}/\code{mc.lbart} function(s). 
For categorical outcomes, use the \code{mbart}/\code{mc.mbart} function(s). 
For survival outcomes, use the \code{surv.bart}/\code{mc.surv.bart} function(s). 
For competing risk outcomes, use the \code{crisk.bart}/\code{mc.crisk.bart} function(s). 
For recurrent event outcomes, use the \code{recur.bart}/\code{mc.recur.bart} function(s). 
For convenience, these functions return the trees (so that you can make
predictions without fitting the model simultaneously) via the \code{predict} function.
}
\author{
Robert McCulloch \email{robert.e.mcculloch@gmail.com}\cr
Maintainer: Rodney Sparapani \email{rsparapa@mcw.edu}
}
%% \references{
%%   Chipman, H., George, E., and McCulloch R. (2010)
%%    Bayesian Additive Regression Trees.
%%    \emph{The Annals of Applied Statistics}, \bold{4,1}, 266-298 <doi:10.1214/09-AOAS285>.

%%    Sparapani, R., Logan, B., McCulloch, R., and Laud, P. (2016)
%%   Nonparametric survival analysis using Bayesian Additive Regression Trees (BART).
%%    \emph{Statistics in Medicine}, \bold{16:2741-53} <doi:10.1002/sim.6893>.
%% }
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{ package }

