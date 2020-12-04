
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017-2018 Robert McCulloch and Rodney Sparapani
## mc.pwbart

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

mc.pwbart = function(
   x.test,		#x matrix to predict at
   treedraws,		#$treedraws from wbart
   mu=0,		#mean to add on
   mc.cores=2L,
   transposed=FALSE,
   dodraws=TRUE,
   nice=19L
)
{
    if(.Platform$OS.type!='unix')
        stop('parallel::mcparallel/mccollect do not exist on windows')

    if(!transposed) x.test <- t(bartModelMatrix(x.test))

    p <- length(treedraws$cutpoints)

    if(p!=nrow(x.test))
        stop(paste0('The number of columns in x.test must be equal to ', p))

    mc.cores.detected <- detectCores()

    if(!is.na(mc.cores.detected) && mc.cores>mc.cores.detected) mc.cores <- mc.cores.detected

    K <- ncol(x.test)
    if(K<mc.cores) mc.cores=K
    k <- K%/%mc.cores-1
    j <- K
    for(i in 1:mc.cores) {
        if(i==mc.cores) h <- 1
        else h <- j-k
        ##print(c(i=i, h=h, j=j))
        parallel::mcparallel({psnice(value=nice);
            pwbart(matrix(x.test[ , h:j], nrow=p, ncol=j-h+1), treedraws, mu, 1, TRUE)},
            silent=(i!=1))
        j <- h-1
    }

    ## K <- ncol(x.test)
    ## k <- ceiling(K/mc.cores)
    ## h <- K-k

    ## parallel::mcparallel({psnice(value=nice);
    ##     pwbart(trees, x.test[ , max(1, h):K], mu, 1, TRUE)})

    ## if(mc.cores>1) for(i in 1:(mc.cores-1)) {
    ##     parallel::mcparallel({psnice(value=nice);
    ##         pwbart(trees, x.test[ , max(1, (h-k)):(h-1)], mu, 1, TRUE)},
    ##         silent=TRUE)
    ##     h <- h-k
    ## }

    pred.list <- parallel::mccollect()
    pred <- pred.list[[1]]
    type=class(pred)[1]
    if(type=='list') pred <- pred[[1]]
    else if(type!='matrix') return(pred.list) ## likely error messages

    if(mc.cores>1) for(i in 2:mc.cores) {
            if(type=='list') pred <- cbind(pred, pred.list[[i]][[1]])
            else pred <- cbind(pred, pred.list[[i]])
        }
    ##if(mc.cores>1) for(i in 2:mc.cores) pred <- cbind(pred, pred.list[[i]])

    if(dodraws) return(pred)
    else return(apply(pred, 2, mean))
    ## if(dodraws) return(pred+mu)
    ## else return(apply(pred, 2, mean)+mu)
}
