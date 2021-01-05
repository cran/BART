\name{leukemia}
\docType{data}
\alias{leukemia}
\title{Bone marrow transplantation for leukemia and multi-state models}

\description{137 patients with acute myelocytic leukemia (AML) and acute
  lymphoblastic leukemia (ALL) were given oral busulfan (Bu) 4 mg/kg on
  each of 4 days and intravenous cyclophosphamide (Cy) 60 mg/kg on each
  of 2 days (BuCy2) followed by allogeneic bone marrow transplantation
  from an HLA-identical or one antigen disparate sibling.}

\usage{data(leukemia)}
\format{A data frame with 137 subjects on the following 22 variables.
\describe{
\item{\code{G}}{Disease Group (1=ALL, 2=AML Low Risk in first
  remission, 3=AML High Risk not in first remission)}
\item{\code{TD}}{Time To Death Or On Study Time}
\item{\code{TB}}{Disease Free Survival Time (Time To Relapse, Death Or End Of Study)}
\item{\code{D}}{Death Indicator (0=Alive, 1=Dead)}
\item{\code{R}}{Relapse Indicator (0=Disease Free, 1=Relapsed)}
\item{\code{B}}{Disease Free Survival Indicator (0=Alive and Disease
  Free, 1=Dead or Relapsed)}
\item{\code{TA}}{Time To Acute Graft-Versus-Host Disease (GVHD)}
\item{\code{A}}{Acute GVHD Indicator (0=Never Developed Acute GVHD, 1=Developed Acute GVHD)}
\item{\code{TC}}{Time To Chronic Graft-Versus-Host Disease (GVHD)}
\item{\code{C}}{Chronic GVHD Indicator (0=Never Developed Chronic GVHD, 1=Developed Chronic GVHD)}
\item{\code{TP}}{Time of Platelets Returning to Normal Levels}
\item{\code{P}}{Platelet Recovery Indicator (0=Platelets Never Returned
  to Normal, 1=Platelets Returned To Normal)}
\item{\code{X1}}{Patient Age In Years}
\item{\code{X2}}{Donor Age In Years}
\item{\code{X3}}{Patient Gender (0=female, 1=male)}
\item{\code{X4}}{Donor Gender (0=female, 1=male)}
\item{\code{X5}}{Patient Cytomegalovirus (CMV) Immune Status (0=CMV Negative, 1=CMV Positive)}
\item{\code{X6}}{Donor Cytomegalovirus (CMV) Immune Status (0=CMV Negative, 1=CMV Positive)}
\item{\code{X7}}{Waiting Time to Transplant In Days}
\item{\code{X8}}{AML Patients with Elevated Risk By
  French-American-British (FAB) Classification (0=Not AML/Elevated, 1=FAB M4 Or M5 with AML)}
\item{\code{X9}}{Hospital (1=The Ohio State University in Columbus,
  2=Alfred in Melbourne, 3=St. Vincent in Sydney, 4=Hahnemann University
in Philadelphia)}
\item{\code{X10}}{Methotrexate Used as a Graft-Versus-Host Disease Prophylactic
  (0=No, 1=Yes)}
}
}

\source{
  Klein J. and Moeschberger M.L. (2003)
  \emph{Survival Analysis: Techniques for Censored and Truncated Data},  
  New York: Springer-Verlag, 2nd Ed., Section 1.3.
  \url{https://www.mcw.edu/-/media/MCW/Departments/Biostatistics/datafromsection13.txt}
 }

\references{
  Copelan E., Biggs J., Thompson J., Crilley P., Szer J., Klein, J., Kapoor N., Avalos, B., Cunningham I. and Atkinson, K. (1991)
  "Treatment for acute myelocytic leukemia with allogeneic bone marrow transplantation following preparation with BuCy2".
  \emph{Blood}, \bold{78(3)}, pp.838-843.
}

\keyword{datasets}
