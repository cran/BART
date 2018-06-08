
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
    x.train = matrix(0,0,0),
    y.train=NULL, times=NULL, delta=NULL,
    x.test = matrix(0,0,0),
    x.test.nogrid = FALSE, ## you may not need the whole grid
    sparse=FALSE, theta=0, omega=1,
    a=0.5, b=1, augment=FALSE, rho=NULL,
    xinfo=matrix(0,0,0), usequants=FALSE,
    ##cont=FALSE,
    rm.const=TRUE, type='pbart',
    ntype=as.integer(
        factor(type, levels=c('wbart', 'pbart', 'lbart'))),
    k = 2, ## BEWARE: do NOT use k for other purposes below
    power = 2, base = 0.95,
    offset = NULL, tau.num=c(NA, 3, 6)[ntype],
    ##binaryOffset = NULL,
    ntree = 50L, numcut = 100L,
    ndpost = 1000L, nskip = 250L,
    keepevery = 10L,
    ##nkeeptrain=ndpost, nkeeptest=ndpost,
    ##nkeeptestmean=ndpost,
    ##nkeeptreedraws=ndpost,
    printevery=100L,
    ##treesaslists=FALSE,
    keeptrainfits=TRUE,
    seed = 99L,    ## only used by mc.recur.bart
    mc.cores = 2L, ## ditto
    nice=19L       ## ditto
    )
{
    if(is.na(ntype) || ntype==1)
        stop("type argument must be set to either 'pbart' or 'lbart'")

    x.train <- bartModelMatrix(x.train)
    x.test <- bartModelMatrix(x.test)

    if(length(rho)==0) rho <- ncol(x.train)

    if(length(y.train)==0) {
        ## if(length(binaryOffset)==0) {
        ##     lambda <- sum(delta, na.rm=TRUE)/
        ##         sum(apply(times, 1, max, na.rm=TRUE))
        ##     binaryOffset <- qnorm(1-exp(-lambda))
        ## }

        recur <- recur.pre.bart(times, delta, x.train)
        ##recur <- recur.pre.bart(times, delta, x.train, x.test)

        y.train <- recur$y.train
        x.train <- recur$tx.train
        x.test  <- recur$tx.test

        times   <- recur$times
        K       <- recur$K
    }
    else {
        times <- unique(sort(x.train[ , 1]))
        K     <- length(times)
    }

    ##if(length(binaryOffset)==0) binaryOffset <- qnorm(mean(y.train))

    ## if(type=='pbart') call <- pbart
    ## else if(type=='lbart') {
    ##     ##binaryOffset <- 0
    ##     call <- lbart
    ## }

    post <- gbart(x.train=x.train, y.train=y.train,
                  x.test=x.test, type=type,
                  sparse=sparse, theta=theta, omega=omega,
                  a=a, b=b, augment=augment, rho=rho,
                  k=k, power=power, base=base,
                  xinfo=xinfo, usequants=usequants,
                  ##cont=cont,
                  rm.const=rm.const,
                  offset=offset, tau.num=tau.num,
                  ##binaryOffset=binaryOffset,
                  ntree=ntree, numcut=numcut,
                  ndpost=ndpost, nskip=nskip,
                  keepevery=keepevery,
                  ##nkeeptrain=nkeeptrain, nkeeptest=nkeeptest,
                  ##nkeeptestmean=nkeeptestmean,
                  ##nkeeptreedraws=nkeeptreedraws,
                  printevery=printevery)

    if(type!=attr(post, 'class')) return(post)

    ##post$binaryOffset <- binaryOffset
    post$times <- times
    post$K <- K
    post$tx.train <- x.train
    post$type <- type

    if(keeptrainfits) {
        ## if(type=='pbart') post$haz.train <- pnorm(post$yhat.train)
        ## else if(type=='lbart') post$haz.train <- plogis(post$yhat.train)

        post$haz.train <- post$prob.train
        post$cum.train <- post$haz.train

        H <- nrow(x.train)

        for(h in 1:H) {
            j <- which(x.train[h, 1]==times) ## for grid points only

            if(j==1) post$haz.train[ , h] <- post$haz.train[ , h]/times[1]
            else {
                post$haz.train[ , h] <- post$haz.train[ , h]/(times[j]-times[j-1])
                post$cum.train[ , h] <- post$cum.train[ , h-1]+post$cum.train[ , h]
            }
        }

        post$haz.train.mean <- apply(post$haz.train, 2, mean)
        post$cum.train.mean <- apply(post$cum.train, 2, mean)
    }

    if(length(x.test)>0) { ## this should always be the case
        post$tx.test <- x.test

        ## if(type=='pbart') post$haz.test <- pnorm(post$yhat.test)
        ## else if(type=='lbart') post$haz.test <- plogis(post$yhat.test)

        if(!x.test.nogrid) {
            post$haz.test <- post$prob.test
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

            post$haz.test.mean <- apply(post$haz.test, 2, mean)
            post$cum.test.mean <- apply(post$cum.test, 2, mean)
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
