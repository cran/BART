
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


recur.bart <- function(
    x.train = matrix(0.0, 0L, 0L),
    y.train=NULL, times=NULL, delta=NULL,
    x.test = matrix(0.0, 0L, 0L),
    x.test.nogrid = FALSE, ## you may not need the whole grid
    k = 2.0, ## BEWARE: do NOT use k for other purposes below
    power = 2.0, base = 0.95,
    binaryOffset = NULL, ##M=1,
    ntree = 50L, numcut = 100L,
    ndpost = 1000L, nskip = 250L,
    keepevery = 10L,
    nkeeptrain=ndpost, nkeeptest=ndpost,
    nkeeptestmean=ndpost, nkeeptreedraws=ndpost,
    printevery=100L,
    treesaslists=FALSE, keeptrainfits=TRUE,
    seed = 99L,    ## only used by mc.recur.bart
    mc.cores = 2L, ## ditto
    nice=19L       ## ditto
    )
{
    if(length(y.train)==0) {
        if(length(binaryOffset)==0) {
            lambda <- sum(delta)/sum(apply(times, 1, max))
            binaryOffset <- qnorm(1-exp(-lambda))
        }

        recur <- recur.pre.bart(times, delta, x.train, x.test)

        y.train <- recur$y.train
        x.train <- recur$tx.train
        x.test  <- recur$tx.test

        times   <- recur$times
        K       <- recur$K
    }
    else {
        if(length(binaryOffset)==0) binaryOffset <- 0

        times <- unique(sort(x.train[ , 1]))
        K     <- length(times)
    }

    post <- pbart(x.train=x.train, y.train=y.train, x.test=x.test,
                  k=k, power=power, base=base,
                  binaryOffset=binaryOffset, ##M=M,
                  ntree=ntree, numcut=numcut,
                  ndpost=ndpost, nskip=nskip,
                  keepevery=keepevery, nkeeptrain=nkeeptrain,
                  nkeeptest=nkeeptest, nkeeptestmean=nkeeptestmean,
                  nkeeptreedraws=nkeeptreedraws, printevery=printevery,
                  treesaslists=treesaslists)

    if(attr(post, 'class')!='pbart') return(post)

    post$binaryOffset <- binaryOffset
    post$times <- times
    post$K <- K
    post$tx.train <- x.train

    ## training grid could be incomplete due to death and/or forced zeros

    ## if(length(x.test)==0) {
    ##     post$cum.train <- pnorm(post$yhat.train)
    ##     post$haz.train <- post$cum.train

    ##     H <- nrow(x.train)

    ##     for(h in 1:H) {
    ##         j <- which(x.train[h, 1]==times) ## for grid points only

    ##         if(j==1) post$haz.train[ , h] <- post$haz.train[ , h]/times[1]
    ##         else {
    ##             post$haz.train[ , h] <- post$haz.train[ , h]/(times[j]-times[j-1])
    ##             post$cum.train[ , h] <- post$cum.train[ , h-1]+post$cum.train[ , h]
    ##         }
    ##     }
    ## }
    ## else {

    if(length(x.test)>0) { ## this should always be the case
        post$tx.test <- x.test

        post$haz.test <- pnorm(post$yhat.test)

        if(!x.test.nogrid) {
            post$cum.test <- post$haz.test

            H <- nrow(x.test)

            for(h in 1:H) {
                j <- which(x.test[h, 1]==times) ## for grid points only

                if(j==1) post$haz.test[ , h] <- post$haz.test[ , h]/times[1]
                else {
                    post$haz.test[ , h] <- post$haz.test[ , h]/(times[j]-times[j-1])
                    post$cum.test[ , h] <- post$cum.test[ , h-1]+post$cum.test[ , h]
                }
            }
        }

        ## H <- nrow(x.test)
        ## time <- 0

        ## for(h in 1:H) {
        ##     prev <- time
        ##     time <- x.test[h, 1]

        ##     if(time==post$times[1]) post$haz.test[ , h] <- post$haz.test[ , h]/time
        ##     else {
        ##         post$haz.test[ , h] <- post$haz.test[ , h]/(time-prev)
        ##         post$cum.test[ , h] <- post$cum.test[ , h-1]+post$cum.test[ , h]
        ##     }
        ## }

    }

    attr(post, 'class') <- 'recurbart'

    return(post)
}
