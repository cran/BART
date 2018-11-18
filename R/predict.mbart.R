
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

predict.mbart <- function(object, newdata, mc.cores=1,
                          openmp=(mc.cores.openmp()>0), ...) {

    ##if(class(newdata) != "matrix") stop("newdata must be a matrix")

    p <- length(object$treedraws$cutpoints)

    if(p!=ncol(newdata))
        stop(paste0('The number of columns in newdata must be equal to ', p))

    if(.Platform$OS.type == "unix") mc.cores.detected <- detectCores()
    else mc.cores.detected <- NA

    if(!is.na(mc.cores.detected) && mc.cores>mc.cores.detected)
        mc.cores <- mc.cores.detected

    if(.Platform$OS.type != "unix" || openmp || mc.cores==1) call <- pwbart
    else call <- mc.pwbart

##call <- predict.gbart
    
    ##return(call(newdata, object$treedraws, mc.cores=mc.cores, mu=object$binaryOffset, ...))

    K <- object$K
    L <- K-1
    pred <- as.list(1:L)
    trees <- object$treedraws$trees
    
    for(h in 1:L) {
        ## eval(parse(text=paste0('object$treedraws$trees=',
        ##                        'object$treedraws$tree', h)))
        object$treedraws$trees <- trees[[h]]
        pred[[h]] <- list(yhat.test=call(newdata, object$treedraws,
                                         mc.cores=mc.cores,
                                         mu=object$offset[h], ...))
        ## predict.gbart testing
        ## pred[[h]] <- call(object, newdata, mc.cores=mc.cores,
        ##                   openmp=openmp, type=object$type)
    }
    H <- dim(pred[[1]]$yhat.test)
    ndpost <- H[1]
    np <- H[2]
    res <- list()
    res$yhat.test <- matrix(nrow=ndpost, ncol=K*np)
    res$prob.test <- matrix(nrow=ndpost, ncol=K*np)
    res$comp.test <- matrix(nrow=ndpost, ncol=K*np)
    
    for(i in 1:np) {
        for(j in 1:K) { 
            h <- (i-1)*K+j
            if(j<K) {
                res$yhat.test[ , h] <- pred[[j]]$yhat.test[ , i]
                if(object$type=='pbart')
                    res$prob.test[ , h] <- pnorm(res$yhat.test[ , h])
                else if(object$type=='lbart')
                    res$prob.test[ , h] <- plogis(res$yhat.test[ , h])
                if(j==1) res$comp.test[ , h] <- 1-res$prob.test[ , h]
                else {
                    res$comp.test[ , h] <- res$comp.test[ , h-1]*
                         (1-res$prob.test[ , h])
                    res$prob.test[ , h] <- res$comp.test[ , h-1]*
                        res$prob.test[ , h] 
                }
            } else res$prob.test[ , h] <- res$comp.test[ , h-1]
        }
    }

    res$prob.test.mean <- apply(res$prob.test, 2, mean)
    res$yhat.test.mean <- NULL
    res$comp.test <- NULL
    res$K <- K
    res$offset <- object$offset
    attr(res, 'class') <- 'mbart'
    
    return(res)
}

