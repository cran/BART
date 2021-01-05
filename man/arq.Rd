\name{arq}
\alias{arq}
\docType{data}
\title{NHANES 2009-2010 Arthritis Questionnaire}

\description{
  This data set was created from the National Health and Nutrition Examination Survey (NHANES)
  2009-2010 Arthritis Questionnaire.  
}

\usage{
data(arq)
}

\details{ We have two outcomes of interest.  Chronic neck pain: Yes
  \code{arq010a=1} vs.\ No \code{arq010a=0}.  Chronic lower-back/buttock
  pain: Yes \code{arq010de=1} vs.\ No \code{arq010de=0}.  \code{seqn} is
  a unique survey respondent identifier.  \code{wtint2yr} is the survey
  sampling weight.  \code{riagendr} is gender: 1 for males, 2 for
  females.  \code{ridageyr} is age in years.  There are several
  anthropometric measurements: \code{bmxwt}, weight in kg; \code{bmxht},
  height in cm; \code{bmxbmi}, body mass index in kg/\eqn{m^2}; and
  \code{bmxwaist}, waist circumference in cm.  The data was subsetted
  to ensure non-missing values of these variables.
}

\references{
  National Health and Nutrition Examination Survey (NHANES)
  2009-2010 Arthritis Questionnaire.
  \url{https://wwwn.cdc.gov/nchs/nhanes/2009-2010/ARQ_F.htm}
}

\keyword{datasets}
