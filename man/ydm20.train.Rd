\name{ydm20.train}
\alias{ydm20.train}
\alias{ydm20.test}
\docType{data}
\title{
A data set used in example of \code{recur.bart}. 
}
\description{
  
Two vectors containing the training and testing set outcomes for a 20\%
random sample for a real data example of recurrent events survival
analysis.  There are 100 patients in the cohort: 50 in the training set
and 50 in the testing set.   See the Reference below (and the References
therein) for more detailed information; a brief synopsis follows.

\code{ydm20.train} contains the training set only.  \code{ydm20.test} is
provided for completeness; it contains both the training set and the
testing set.  There are 798 unique time points so there are 50*798=39900
rows of the training set followed by 50*798=39900 rows of the testing
set.

The full data sets, \code{ydm.train} and \code{ydm.test}, can be
obtained online at
%\url{http://www.mcw.edu/FileLibrary/Groups/Biostatistics/TechReports/TechReports5175/tr064.tar}.
\url{https://www.mcw.edu/-/media/MCW/Departments/Biostatistics/Tech-Report-64-April-2017-data-file.tar?la=en}.
There are 488 patients in the full cohort: 235 in the training set and
253 in the testing set.

\code{ydm.train} contains the training set only.  \code{ydm.test}
contains both the training set and the testing set.  There are 798
unique time points so there are 235*798=187530 rows of the training set
followed by 253*798=201894 rows of the testing set. 

}
\usage{
data(ydm20.train)
data(ydm20.test)
}
\seealso{xdm20.train}
\references{
Sparapani, Rein, Tarima, Jackson, Meurer (2020).  Non-parametric recurrent
events analysis with BART and an application to the hospital admissions
of patients with diabetes. \emph{Biostatistics} doi:10.1093/biostatistics/kxy032
%% Sparapani R., Rein L., Tarima S., Jackson T., Meurer J. (2017)
%% Nonparametric recurrent events analysis with BART
%% and an application to the hospital admissions of patients with diabetes.
%% \emph{MCW Biostatistics Technical Report 64}
%% \url{http://www.mcw.edu/FileLibrary/Groups/Biostatistics/TechReports/TechReports5175/tr064.pdf}.
}
\examples{
data(ydm20.train)
data(ydm20.test)
table(ydm20.train)
table(ydm20.test)
}
\keyword{datasets}
