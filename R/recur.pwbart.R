
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

recur.pwbart <- function(
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

    if(class(pred$yhat.test)!='matrix') return(pred$yhat.test)

    x.test <- t(x.test)
    pred$tx.test <- x.test
    times <- unique(sort(x.test[ , 1]))
    pred$times <- times
    K <- length(times)
    pred$K <- K

    if(type=='pbart') pred$prob.test <- pnorm(pred$yhat.test)
    else if(type=='lbart') pred$prob.test <- plogis(pred$yhat.test)

    pred$haz.test <- pred$prob.test
    pred$cum.test <- pred$haz.test

    H <- nrow(x.test)

    for(h in 1:H) {
        j <- which(x.test[h, 1]==times) ## for grid points only

        if(j==1) pred$haz.test[ , h] <- pred$haz.test[ , h]/times[1]
        else {
            pred$haz.test[ , h] <- pred$haz.test[ , h]/(times[j]-times[j-1])
            pred$cum.test[ , h] <- pred$cum.test[ , h-1]+pred$cum.test[ , h]
        }
    }

    pred$prob.test.mean <- apply(pred$prob.test, 2, mean)
    pred$haz.test.mean <- apply(pred$haz.test, 2, mean)
    pred$cum.test.mean <- apply(pred$cum.test, 2, mean)
    
    pred$binaryOffset <- binaryOffset
    attr(pred, 'class') <- 'recurbart'

    return(pred)
}
