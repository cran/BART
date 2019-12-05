
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

surv.pwbart <- function(
   x.test,		#x matrix to predict at with time points expanded
   treedraws,		#$treedraws for from surv.bart/mc.surv.bart
   binaryOffset=0,	#mean to add on
   mc.cores=1L,
   type='pbart',
   transposed=FALSE,
   nice=19L             # mc.surv.pwbart only
)
{
    if(!transposed) x.test <- t(bartModelMatrix(x.test))

    p <- length(treedraws$cutpoints)

    if(p!=nrow(x.test))
        stop(paste0('The number of columns in x.test must be equal to ', p))

    pred <- list()

    pred$yhat.test <- pwbart(x.test, treedraws, binaryOffset, mc.cores, TRUE)

    if(class(pred$yhat.test)[1]!='matrix') return(pred$yhat.test)

    x.test <- t(x.test)
    pred$tx.test <- x.test
    times <- unique(sort(x.test[ , 1]))
    pred$times <- times
    K <- length(times)
    pred$K <- K

    H <- nrow(x.test)/K ## the number of different settings

    if(type=='pbart') pred$prob.test <- pnorm(pred$yhat.test)
    else if(type=='lbart') pred$prob.test <- plogis(pred$yhat.test)

    pred$surv.test <- 1-pred$prob.test

    for(h in 1:H)
        for(j in 2:K) {
                      l <- K*(h-1)+j

                      pred$surv.test[ , l] <- pred$surv.test[ , l-1]*pred$surv.test[ , l]
        }

    pred$surv.test.mean <- apply(pred$surv.test, 2, mean)

    pred$binaryOffset <- binaryOffset
    attr(pred, 'class') <- 'survbart'

    return(pred)
}
