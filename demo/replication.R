
library(BART)

## replication script for JSS article

## directory to create graphics files in, if any
## if option not specified, then default assumed to be 'NONE'
options(figures='.')
##options(figures='../vignettes/figures')

## for single-threading, specify one core
## Windows lacks forking (and generally lacks OpenMP)
## so single-threading only
## for multi-threading, specify the number of cores
## a Unix-like OS provides forking for multi-threading
## (and often OpenMP is available as well)

if(.Platform$OS.type=='unix') {
    ## there are diminishing returns so often 8 cores is sufficient
    options(mc.cores=min(8, parallel::detectCores()))
} else {
    options(mc.cores=1)
}

## Section 3, The Boston Housing Data including Figures 1-6
source(system.file('demo/boston.R', package='BART'))

## Section 4.2, Probit BART Example: Chronic Pain and Obesity
## Figure 7
source(system.file('demo/nhanes.pbart1.R', package='BART'))
## Figure 8
source(system.file('demo/nhanes.pbart2.R', package='BART'))

## Section 4.4, Multinomial BART Example: Alligator Food Preference
## Figure 9
source(system.file('demo/alligator.R', package='BART'))

## Section 4.5, Convergence Diagnostics for Binary and Categorical Outcomes
## Figures 10-12
source(system.file('demo/geweke.pbart2.R', package='BART'))

## Section 4.6, BART and Variable Selection
## Figure 13
source(system.file('demo/sparse.pbart.R', package='BART'))

## Section 5.1, Survival Analysis with BART Example: Advanced Lung Cancer
## Figure 14
source(system.file('demo/lung.surv.bart.R', package='BART'))

## Section 5.3, Competing Risks with BART Example: Liver Transplants
## Figure 15
source(system.file('demo/liver.crisk.bart.R', package='BART'))

## Section 5.4, Recurrent Events with BART Example: Bladder Tumors}
## Figures 16-18
source(system.file('demo/bladder.recur.bart.R', package='BART'))
