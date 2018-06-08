
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017 Robert McCulloch and Rodney Sparapani

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## https://www.R-project.org/Licenses/GPL-2


## you call this function before crisk.bart()
## this function takes traditional time/delta
## competing risk variables and regressors (if any)
## and it constructs the corresponding
## tx.train, y.train, y.train2, and tx.test appropriate
## for use with pbart()

crisk.pre.bart <- function(
                      times,
                      ## vector of survival times

                      delta,
                      ## vector of event indicators
                      ## 0=censoring, 1=cause 1, 2=cause 2

                      x.train=NULL,
                      ## matrix of covariate regressors for cause 1
                      ## can be NULL, i.e. KM analog

                      x.test=NULL,
                      ## matrix of covariate regressors for cause 1

                      x.train2=x.train,
                      ## matrix of covariate regressors for cause 2
                      ## can be NULL, i.e. KM analog

                      x.test2=x.test,
                      ## matrix of covariate regressors for cause 2

                      K=NULL
                      ## if specified, then use K quantiles for time grid
                      ) {
    ## currently does not handle time dependent Xs
    ## can be extended later
    ## most likely via the alternative counting process notation

    check <- unique(sort(delta))

    C <- length(check)

    ## with gap times, there is no censoring per se, only 1s and 2s
    if(!(C %in% 2:3) || (C==3 && !all(check==0:2)) || (C==2 && !all(check==1:2)))
        stop('delta must be coded as: 0(censored), 1(cause 1) or 2(cause 2)')

    if(length(x.train)>0 && length(x.train2)>0 && nrow(x.train)!=nrow(x.train2))
        stop('number of rows in x.train and x.train2 must be equal')

    if(length(x.test)>0 && length(x.test2)>0 && nrow(x.test)!=nrow(x.test2))
        stop('number of rows in x.test and x.test2 must be equal')

    pre <- surv.pre.bart(times=times, 1*(delta==1), K=K,
                         x.train=x.train, x.test=x.test)

    pre$cond <- which(pre$y.train==0)

    pre2 <- surv.pre.bart(times=times, 1*(delta==2), K=K,
                         x.train=x.train2, x.test=x.test2)

    pre$tx.train2 <- pre2$tx.train
    pre$tx.test2 <- pre2$tx.test
    pre$y.train2 <- pre2$y.train
    ##pre$binaryOffset2 <- pre2$binaryOffset

    return(pre)
}
