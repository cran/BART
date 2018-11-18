
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2018 Robert McCulloch and Rodney Sparapani

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


crisk2.bart <- function(
    x.train=matrix(0,0,0), y.train=NULL, x.train2=x.train,
    y.train2=NULL, times=NULL, delta=NULL, K=NULL,
    x.test=matrix(0,0,0), x.test2=x.test,
    sparse=FALSE, theta=0, omega=1,
    a=0.5, b=1, augment=FALSE, rho=NULL, rho2=NULL,
    xinfo=matrix(0,0,0), xinfo2=matrix(0,0,0), usequants=FALSE,
    ##cont=FALSE,
    rm.const=TRUE, type='pbart',
    ntype=as.integer(
        factor(type, levels=c('wbart', 'pbart', 'lbart'))),
    k = 2, ## BEWARE: do NOT use k for other purposes below
    power = 2, base = 0.95,
    offset = NULL, offset2 = NULL,
    tau.num=c(NA, 3, 6)[ntype],
    ##binaryOffset = NULL, binaryOffset2 = NULL,
    ntree = 50L, numcut = 100L,
    ndpost = 1000L, nskip = 250L,
    keepevery = 10L,
    ##nkeeptrain=ndpost, nkeeptest=ndpost, nkeeptreedraws=ndpost,
    printevery=100L,
    ##keeptrainfits=TRUE,
    id = NULL,
    seed=99,    ## mc.crisk2.bart only
    mc.cores=2, ## mc.crisk2.bart only
    nice=19L    ## mc.crisk2.bart only
)
{
    if(is.na(ntype) || ntype==1)
        stop("type argument must be set to either 'pbart' or 'lbart'")

    x.train2 <- bartModelMatrix(x.train2)
    x.test2 <- bartModelMatrix(x.test2)
    x.train <- bartModelMatrix(x.train)
    x.test <- bartModelMatrix(x.test)

    if(length(rho)==0) rho=ncol(x.train)
    if(length(rho2)==0) rho2=ncol(x.train2)

    if(length(y.train)==0) {
        pre <- surv.pre.bart(times, delta, x.train, x.test, K=K)
        pre2 <- surv.pre.bart(times, delta, x.train2, x.test2, K=K)

        y.train <- pre$y.train
        x.train <- pre$tx.train
        x.test  <- pre$tx.test
        y.train2 <- 1*(y.train[y.train>0]==1)
        x.train2 <- cbind(pre2$tx.train[y.train>0, ])
        x.test2  <- pre2$tx.test
        y.train <- 1*(y.train>0)

        times   <- pre$times
        K       <- pre$K
    }
    else {
        ## if(length(x.train)>0 & length(x.train2)>0 &
        ##    nrow(x.train)!=nrow(x.train2))
        ##     stop('number of rows in x.train and x.train2 must be equal')

        if(length(x.test)>0 & length(x.test2)>0 & nrow(x.test)!=nrow(x.test2))
            stop('number of rows in x.test and x.test2 must be equal')

        times <- unique(sort(x.train[ , 1]))
        K     <- length(times)
    }

    if(length(xinfo)==0) {
        temp = bartModelMatrix(x.train2, numcut, usequants=usequants,
                               ##cont=cont,
                               xinfo=xinfo, rm.const=rm.const)
        x.train2 = t(temp$X)
        numcut2 = temp$numcut
        xinfo2 = temp$xinfo
        ## if(length(x.test2)>0)
        ##     x.test2 = t(bartModelMatrix(x.test2[ , temp$rm.const]))
        if(length(x.test2)>0) {
            x.test2 = bartModelMatrix(x.test2)
            x.test2 = t(x.test2[ , temp$rm.const])
        }
        rm.const2 <- temp$rm.const
        rm(temp)

        temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                               ##cont=cont,
                               xinfo=xinfo, rm.const=rm.const)
        x.train = t(temp$X)
        numcut = temp$numcut
        xinfo = temp$xinfo
        ## if(length(x.test)>0)
        ##     x.test = t(bartModelMatrix(x.test[ , temp$rm.const]))
        if(length(x.test)>0) {
            x.test = bartModelMatrix(x.test)
            x.test = t(x.test[ , temp$rm.const])
        }
        rm.const <- temp$rm.const
        rm(temp)

        xinfo2[1, ] <- xinfo[1, ] ## same time grid
        transposed <- TRUE
    }
    else {
        rm.const <- 1:ncol(x.train)
        rm.const2 <- 1:ncol(x.train2)
        transposed <- FALSE
    }

    ## if(type=='pbart') call <- pbart
    ## else if(type=='lbart') {
    ##     ##binaryOffset <- 0
    ##     ##binaryOffset2 <- 0
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
                  ndpost=ndpost, nskip=nskip,
                  keepevery=keepevery, ##nkeeptrain=0,
                  ##nkeeptest=nkeeptest, nkeeptreedraws=nkeeptreedraws,
                  printevery=printevery, transposed=transposed)

    if(type!=attr(post, 'class')) return(post)

    post2 <- gbart(x.train=x.train2, y.train=y.train2,
                   x.test=x.test2, type=type,
                   sparse=sparse, theta=theta, omega=omega,
                   a=a, b=b, augment=augment, rho=rho2,
                   xinfo=xinfo2, usequants=usequants,
                   ##cont=cont,
                   rm.const=rm.const,
                   k=k, power=power, base=base,
                   offset=offset2, tau.num=tau.num,
                   ##binaryOffset=binaryOffset2,
                   ntree=ntree, numcut=numcut2,
                   ndpost=ndpost, nskip=nskip,
                   keepevery=keepevery, ##nkeeptrain=0,
                   ##nkeeptest=nkeeptest,
                   ##nkeeptreedraws=nkeeptreedraws,
                   printevery=printevery, transposed=transposed)

    if(type!=attr(post2, 'class')) return(post2)

    ##post$binaryOffset <- binaryOffset
    post$offset2 <- post2$offset
    ##post$binaryOffset2 <- post2$binaryOffset
    post$id <- id
    post$times <- times
    post$K <- K

    if(!transposed) {
        post$tx.train <- x.train
        post$tx.train2 <- x.train2
    } else {
        post$tx.train <- t(x.train)
        post$tx.train2 <- t(x.train2)
    }

    post$type <- type
    post$treedraws2 <- post2$treedraws
    post$varcount2 <- post2$varcount
    post$varcount2.mean <- post2$varcount.mean
    post$varprob2 <- post2$varprob
    post$varprob2.mean <- post2$varprob.mean
    post$rm.const <- rm.const
    post$rm.const2 <- rm.const2
    post$yhat.train <- NULL
    post$yhat.train.mean <- NULL

    if(length(x.test)>0) {
        if(!transposed) {
            post$tx.test <- x.test
            post$tx.test2 <- x.test2
        } else {
            post$tx.test <- t(x.test)
            post$tx.test2 <- t(x.test2)
        }

        H <- nrow(post$tx.test)/K ## the number of different settings

        post$yhat.test2 <- post2$yhat.test
        post$prob.test2 <- post2$prob.test

        post$surv.test <- 1-post$prob.test
        post$cif.test <- post$prob.test*post$prob.test2
        post$cif.test2 <- post$prob.test*(1-post$prob.test2)

        for(h in 1:H)
            for(j in 2:K) {
                l <- K*(h-1)+j

                post$cif.test[ , l] <- post$cif.test[ , l-1]+
                    post$surv.test[ , l-1]*post$cif.test[ , l]
                post$cif.test2[ , l] <- post$cif.test2[ , l-1]+
                    post$surv.test[ , l-1]*post$cif.test2[ , l]
                post$surv.test[ , l] <- post$surv.test[ , l-1]*
                    post$surv.test[ , l]
            }

        post$cif.test.mean <- apply(post$cif.test, 2, mean)
        post$cif.test2.mean <- apply(post$cif.test2, 2, mean)
        post$surv.test.mean <- apply(post$surv.test, 2, mean)
    }

    attr(post, 'class') <- 'crisk2bart'

    return(post)
}
