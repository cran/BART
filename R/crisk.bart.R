
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


crisk.bart <- function(
    x.train=matrix(0.0, 0L, 0L), y.train=NULL,
    x.train2=x.train, y.train2=NULL,
    times=NULL, delta=NULL,
    x.test=matrix(0.0, 0L, 0L), x.test2=x.test,
    cond=NULL, k = 2.0, ## BEWARE: do NOT use k for other purposes below
    power = 2.0, base = 0.95,
    binaryOffset = NULL,
    binaryOffset2 = NULL,
    ntree = 50L, numcut = 100L,
    ndpost = 10000L, nskip = 250L,
    keepevery = 10L, 
    nkeeptrain=ndpost%/%keepevery, nkeeptest=ndpost%/%keepevery,
    nkeeptestmean=ndpost%/%keepevery, nkeeptreedraws=ndpost%/%keepevery,
    printevery=100L, 
    treesaslists=FALSE, keeptrainfits=TRUE,
    id = NULL,
    seed=99,    ## mc.crisk.bart only
    mc.cores=2, ## mc.crisk.bart only
    nice=19L    ## mc.crisk.bart only  
)
{
    if(length(y.train)==0) {
        pre <- crisk.pre.bart(times, delta, x.train, x.test, x.train2, x.test2)

        y.train <- pre$y.train
        x.train <- pre$tx.train
        x.test  <- pre$tx.test
        y.train2 <- pre$y.train2
        x.train2 <- pre$tx.train2
        x.test2  <- pre$tx.test2

        times   <- pre$times
        K       <- pre$K

        if(length(cond)==0) cond <- pre$cond
        if(length(binaryOffset)==0) binaryOffset <- pre$binaryOffset
        if(length(binaryOffset2)==0) binaryOffset2 <- pre$binaryOffset2
    }
    else {
        if(length(x.train)>0 & length(x.train2)>0 & nrow(x.train)!=nrow(x.train2))
            stop('number of rows in x.train and x.train2 must be equal')

        if(length(x.test)>0 & length(x.test2)>0 & nrow(x.test)!=nrow(x.test2))
            stop('number of rows in x.test and x.test2 must be equal')

        if(length(binaryOffset)==0) binaryOffset <- 0
        if(length(binaryOffset2)==0) binaryOffset2 <- 0

        times <- unique(sort(x.train[ , 1]))
        K     <- length(times)
    }

    post <- pbart(x.train=x.train, y.train=y.train, x.test=x.test,
                  k=k, power=power, base=base,
                  binaryOffset=binaryOffset,
                  ntree=ntree, numcut=numcut,
                  ndpost=ndpost, nskip=nskip,
                  keepevery=keepevery, nkeeptrain=0*nkeeptrain,
                  nkeeptest=nkeeptest, nkeeptestmean=nkeeptestmean,
                  nkeeptreedraws=nkeeptreedraws, printevery=printevery,
                  treesaslists=treesaslists)

    post2 <- pbart(x.train=as.matrix(x.train2[cond, ]),
                   y.train=y.train2[cond], x.test=x.test2,
                   k=k, power=power, base=base,
                   binaryOffset=binaryOffset2,
                   ntree=ntree, numcut=numcut,
                   ndpost=ndpost, nskip=nskip,
                   keepevery=keepevery, nkeeptrain=0*nkeeptrain,
                   nkeeptest=nkeeptest, nkeeptestmean=nkeeptestmean,
                   nkeeptreedraws=nkeeptreedraws, printevery=printevery,
                   treesaslists=treesaslists)

    post$binaryOffset <- binaryOffset
    post$binaryOffset2 <- binaryOffset2
    post$id <- id
    post$times <- times
    post$K <- K
    post$tx.train <- x.train
    post$tx.train2 <- x.train2
    post$cond <- cond
    post$treedraws2 <- post2$treedraws
    post$varcount2 <- post2$varcount
    post$yhat.train <- NULL
    post$yhat.train.mean <- NULL
    
    if(length(x.test)>0) {
        post$tx.test <- x.test
        post$tx.test2 <- x.test2
        H <- nrow(x.test)/K ## the number of different settings

        post$yhat.test2 <- post2$yhat.test

        post$prob.test <- pnorm(post$yhat.test)
        post$prob.test2 <- pnorm(post$yhat.test2)
        post$surv.test <- (1-post$prob.test)*(1-post$prob.test2)
        post$prob.test2 <- (1-post$prob.test)*post$prob.test2
        post$cif.test <- post$prob.test
        post$cif.test2 <- post$prob.test2

        for(h in 1:H) for(j in 2:K) {
                l <- K*(h-1)+j
                
                post$cif.test[ , l] <- post$cif.test[ , l-1]+post$surv.test[ , l-1]*post$cif.test[ , l]
                post$cif.test2[ , l] <- post$cif.test2[ , l-1]+post$surv.test[ , l-1]*post$cif.test2[ , l]
                post$surv.test[ , l] <- post$surv.test[ , l-1]*post$surv.test[ , l]
                      }

        post$cif.test.mean <- apply(post$cif.test, 2, mean)
        post$cif.test2.mean <- apply(post$cif.test2, 2, mean)
        post$surv.test.mean <- apply(post$surv.test, 2, mean)
    }

    return(post)
}
