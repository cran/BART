
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

    ##return(call(newdata, object$treedraws, mc.cores=mc.cores, mu=object$binaryOffset, ...))

    C <- object$C

    pred <- as.list(1:C)

    for(h in 1:C) {
        eval(parse(text=paste0('object$treedraws$trees=',
                               'object$treedraws$tree', h)))

        pred[[h]] <- list(yhat.test=call(newdata, object$treedraws,
                                    mc.cores=mc.cores,
                                    mu=object$binaryOffset[h], ...))
    }

    H <- dim(pred[[1]]$yhat.test)
    np <- H[2]
    res <- list()
    res$yhat.test <- matrix(nrow=H[1], ncol=C*np)
    res$prob.test <- matrix(nrow=H[1], ncol=C*np)
    
    for(i in 1:np) {
        h <- (i-1)*C
        for(j in 1:C) { 
            res$yhat.test[ , h+j] <- pred[[j]]$yhat.test[ , i]
            res$prob.test[ , h+j] <- pnorm(res$yhat.test[ , h+j])
        }
        total <- apply(res$prob.test[ , h+1:C], 1, sum)
        for(j in 1:C) 
            res$prob.test[ , h+j] <- res$prob.test[ , h+j]/total
    }

    res$prob.test.mean <- apply(res$prob.test, 2, mean)
    res$yhat.test.mean <- NULL
    res$C <- C
    res$binaryOffset <- object$binaryOffset
    attr(res, 'class') <- 'mbart'
    
    return(res)
}

