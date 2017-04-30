\name{ACTG175}
\docType{data}
\alias{ACTG175}
\title{AIDS Clinical Trials Group Study 175}
\description{ACTG 175 was a randomized clinical trial to compare monotherapy with zidovudine or didanosine with 
combination therapy with zidovudine and didanosine or zidovudine and zalcitabine in adults infected with the human 
immunodeficiency virus type I whose CD4 T cell counts were between 200 and 500 per cubic millimeter.}
\usage{data(ACTG175)}
\format{A data frame with 2139 observations on the following 27 variables:
\describe{
\item{\code{pidnum}}{patient's ID number}
\item{\code{age}}{age in years at baseline}
\item{\code{wtkg}}{weight in kg at baseline}
\item{\code{hemo}}{hemophilia (0=no, 1=yes)}
\item{\code{homo}}{homosexual activity (0=no, 1=yes)}
\item{\code{drugs}}{history of intravenous drug use (0=no, 1=yes)}
\item{\code{karnof}}{Karnofsky score (on a scale of 0-100)}
\item{\code{oprior}}{non-zidovudine antiretroviral therapy prior to initiation of study treatment (0=no, 1=yes)}
\item{\code{z30}}{zidovudine use in the 30 days prior to treatment initiation (0=no, 1=yes)}
\item{\code{zprior}}{zidovudine use prior to treatment initiation (0=no, 1=yes)}
\item{\code{preanti}}{number of days of previously received antiretroviral therapy}
\item{\code{race}}{race (0=white, 1=non-white)}
\item{\code{gender}}{gender (0=female, 1=male)}
\item{\code{str2}}{antiretroviral history (0=naive, 1=experienced)}
\item{\code{strat}}{antiretroviral history stratification (1='antiretroviral naive', 2='> 1 but \eqn{\leq} 52 weeks of prior
antiretroviral therapy', 3='> 52 weeks')}
\item{\code{symptom}}{symptomatic indicator (0=asymptomatic, 1=symptomatic)}
\item{\code{treat}}{treatment indicator (0=zidovudine only, 1=other therapies)}
\item{\code{offtrt}}{indicator of off-treatment before 96\eqn{\pm}5 weeks (0=no,1=yes)}
\item{\code{cd40}}{CD4 T cell count at baseline}
\item{\code{cd420}}{CD4 T cell count at 20\eqn{\pm}5 weeks}
\item{\code{cd496}}{CD4 T cell count at 96\eqn{\pm}5 weeks (=\code{NA} if missing)}
\item{\code{r}}{missing CD4 T cell count at 96\eqn{\pm}5 weeks (0=missing, 1=observed)}
\item{\code{cd80}}{CD8 T cell count at baseline}
\item{\code{cd820}}{CD8 T cell count at 20\eqn{\pm}5 weeks}
\item{\code{cens}}{indicator of observing the event in \code{days}}
\item{\code{days}}{number of days until the first occurrence of: (i) a decline in CD4 T cell count of at least 50%,
 (ii) an event indicating progression to AIDS, or (iii) death.}
\item{\code{arms}}{treatment arm (0=zidovudine, 1=zidovudine and didanosine, 2=zidovudine and zalcitabine, 
3=didanosine).}
}
}
\details{The variable \code{days} contains right-censored time-to-event observations. The data set includes 
the following post-randomization covariates: CD4 and CD8 T cell count at 20\eqn{\pm}5 
weeks and the indicator of whether or not the patient was taken off-treatment before 96\eqn{\pm}5 weeks.}
\references{Hammer SM, et al. (1996), "A trial comparing nucleoside monotherapy with combination therapy in 
HIV-infected adults with CD4 cell counts from 200 to 500 per cubic millimeter.", \emph{New England Journal of Medicine}, 
335:1081--1090.}
\keyword{datasets}
