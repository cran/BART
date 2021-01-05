\name{bladder}
\docType{data}
\alias{bladder}
\alias{bladder1}
\alias{bladder2}
\title{Bladder Cancer Recurrences}
%\usage{
%bladder
%bladder1
%bladder2
%}
\description{
  
This interesting example is from a clinical trial conducted by the
Veterans Administration Cooperative Urological Research Group.  This
data on recurrence of bladder cancer has been used by many to
demonstrate methodology for recurrent events modelling.  In this study,
all patients had superficial bladder tumors when they entered the
trial. These tumors were removed transurethrally and patients were
randomly assigned to one of three treatments: placebo, thiotepa or
pyridoxine (vitamin B6). Many patients had multiple recurrences of
tumors during the study and new tumors were removed at each visit.  For
each patient, their recurrence time, if any, was measured from the
beginning of treatment.

bladder is the data set that appears most commonly in the literature.
It uses only the 85 subjects with nonzero follow-up who were assigned to
either thiotepa or placebo and only the first four recurrences for any
patient.  The status variable is 1 for recurrence and 0 for everything
else (including death for any reason).  The data set is laid out in the
competing risks format of the paper by Wei, Lin, and Weissfeld (WLW).

bladder1 is the full data set from the study. It contains all three
treatment arms and all recurrences for 118 subjects; the maximum
observed number of recurrences is 9.

bladder2 uses the same subset of subjects as bladder, but formated in
the (start, stop] or Anderson-Gill (AG) style.  Note that in
transforming from the WLW to the AG style data set there is a quite
common programming mistake that leads to extra follow-up time for 12
subjects: all those with follow-up beyond their fourth recurrence.  Over
this extended time these subjects are by definition not at risk for
another event in the WLW data set.

}
\format{
  bladder
  \tabular{ll}{
    id:\tab Patient id\cr
    rx:\tab Treatment 1=placebo  2=thiotepa\cr
    number:\tab Initial number of tumours (8=8 or more)\cr
    size:\tab size (cm) of largest initial tumour\cr
    stop:\tab recurrence or censoring time\cr
    enum:\tab which recurrence (up to 4)\cr
  }
  bladder1
  \tabular{ll}{
    id:\tab Patient id\cr
    treatment:\tab Placebo, pyridoxine (vitamin B6), or thiotepa\cr
    number:\tab Initial number of tumours (8=8 or more)\cr
    size:\tab Size (cm) of largest initial tumour\cr
    recur:\tab Number of recurrences \cr
    start,stop:\tab The start and end time of each time interval\cr
    status:\tab End of interval code, 0=censored, 1=recurrence, \cr
           \tab 2=death from bladder disease, 3=death other/unknown cause\cr
    rtumor:\tab Number of tumors found at the time of a recurrence\cr
    rsize:\tab Size of largest tumor at a recurrence\cr
    enum:\tab Event number (observation number within patient)\cr
  }
  bladder2 
  \tabular{ll}{
    id:\tab Patient id\cr
    rx:\tab Treatment 1=placebo  2=thiotepa\cr
    number:\tab Initial number of tumours (8=8 or more)\cr
    size:\tab size (cm) of largest initial tumour\cr
    start:\tab start of interval (0 or previous recurrence time)\cr
    stop:\tab recurrence or censoring time\cr
    enum:\tab which recurrence (up to 4)\cr
  }
}

\references{
  Byar, DP (1980),
  "The Veterans Administration Study of Chemoprophylaxis for Recurrent Stage I Bladder Tumors: Comparisons of Placebo, Pyridoxine, and Topical Thiotepa,"
  in \emph{Bladder Tumors and Other Topics in Urological Oncology}, eds.
  M Pavone-Macaluso, PH Smith, and F Edsmyn, New York: Plenum, pp. 363-370.

  Andrews DF, Hertzberg AM (1985), 
  DATA: A Collection of Problems from Many Fields for the Student and Research Worker, New York: Springer-Verlag.

  LJ Wei, DY Lin, L Weissfeld (1989),
  Regression analysis of multivariate incomplete failure time data by
  modeling marginal distributions.
  \emph{Journal of the American Statistical Association},
  \bold{84}.
}

\examples{
  data(bladder)
}
\keyword{datasets}
\keyword{survival}
