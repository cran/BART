
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017-2018 Robert McCulloch and Rodney Sparapani

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


surv.bart <- function(
    x.train = matrix(0,0,0),
    y.train=NULL, times=NULL, delta=NULL,
    x.test = matrix(0,0,0),
    K=NULL, events=NULL, ztimes=NULL, zdelta=NULL,
    sparse=FALSE, theta=0, omega=1,
    a=0.5, b=1, augment=FALSE, rho=NULL,
    xinfo=matrix(0,0,0), usequants=FALSE,
    ## cont=FALSE,
    rm.const=TRUE, type='pbart',
    ntype=as.integer(
        factor(type, levels=c('wbart', 'pbart', 'lbart'))),
    k = 2, ## BEWARE: do NOT use k for other purposes below
    power = 2, base = 0.95,
    offset = NULL, tau.num=c(NA, 3, 6)[ntype],
    ##binaryOffset = NULL,
    ntree = 50L, numcut = 100L,
    ndpost = 1000L, nskip = 250L, keepevery = 10L,
    ##nkeeptrain=ndpost, nkeeptest=ndpost,
    ##nkeeptreedraws=ndpost,
    printevery=100L,
    ##keeptrainfits=TRUE,
    id = NULL,     ## only used by surv.bart
    seed = 99L,    ## only used by mc.surv.bart
    mc.cores = 2L, ## ditto
    nice=19L       ## ditto
)
{
    if(is.na(ntype) || ntype==1)
        stop("type argument must be set to either 'pbart' or 'lbart'")

    x.train <- bartModelMatrix(x.train)
    x.test <- bartModelMatrix(x.test)

    if(length(rho)==0) rho=ncol(x.train)

    if(length(y.train)==0) {
        pre <- surv.pre.bart(times, delta, x.train, x.test, K=K,
                             events=events, ztimes=ztimes, zdelta=zdelta)

        y.train <- pre$y.train
        x.train <- pre$tx.train
        x.test  <- pre$tx.test

        times   <- pre$times
        K       <- pre$K

        ##if(length(binaryOffset)==0) binaryOffset <- pre$binaryOffset
    }
    else {
        if(length(unique(sort(y.train)))>2)
            stop('y.train has >2 values; make sure you specify times=times & delta=delta')

        ##if(length(binaryOffset)==0) binaryOffset <- 0

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
                  xinfo=xinfo, usequants=usequants,
                  ##cont=cont,
                  rm.const=rm.const,
                  k=k, power=power, base=base,
                  offset=offset, tau.num=tau.num,
                  ##binaryOffset=binaryOffset,
                  ntree=ntree, numcut=numcut,
                  ndpost=ndpost, nskip=nskip, keepevery=keepevery,
                  ##nkeeptrain=nkeeptrain, nkeeptest=nkeeptest,
                  ##nkeeptestmean=nkeeptestmean,
                  ##nkeeptreedraws=nkeeptreedraws,
                  printevery=printevery)

    if(type!=attr(post, 'class')) return(post)

    ##post$binaryOffset <- binaryOffset
    post$id <- id
    post$times <- times
    post$K <- K
    post$tx.train <- x.train
    post$type <- type

    ## if(keeptrainfits) {
    ##      post$surv.train <- 1-post$prob.train

    ##      H <- nrow(x.train)/K ## the number of different settings

    ##      for(h in 1:H)
    ##          for(j in 2:K) {
    ##              l <- K*(h-1)+j

    ##              post$surv.train[ , l] <-
    ## post$surv.train[ , l-1]*post$surv.train[ , l]
    ##          }

    ##      post$surv.train.mean <- apply(post$surv.train, 2, mean)
    ## }

    if(length(x.test)>0) {
        post$tx.test <- x.test
        H <- nrow(x.test)/K ## the number of different settings

        post$surv.test <- 1-post$prob.test
        ## if(type=='pbart') post$surv.test <- 1-pnorm(post$yhat.test)
        ## else if(type=='lbart') post$surv.test <- 1-plogis(post$yhat.test)

        for(h in 1:H)
            for(j in 2:K) {
                l <- K*(h-1)+j

                post$surv.test[ , l] <-
                    post$surv.test[ , l-1]*post$surv.test[ , l]
            }

        post$surv.test.mean <- apply(post$surv.test, 2, mean)
    }

    attr(post, 'class') <- 'survbart'

    return(post)
}
