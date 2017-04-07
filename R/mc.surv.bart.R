
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

## run BART and generate survival in parallel

mc.surv.bart <- function(
    x.train = matrix(0.0, 0L, 0L),
    y.train=NULL, times=NULL, delta=NULL,
    x.test = matrix(0.0, 0L, 0L),
    k = 2.0, ## BEWARE: do NOT use k for other purposes below
    power = 2.0, base = 0.95,
    binaryOffset = NULL,
    ntree = 50L, numcut = 100L,
    ndpost = 10000L, nskip = 250L,
    keepevery = 10L, 
    nkeeptrain=ndpost%/%keepevery, nkeeptest=ndpost%/%keepevery,
    nkeeptestmean=ndpost%/%keepevery, nkeeptreedraws=ndpost%/%keepevery,
    printevery=100L, 
    treesaslists=FALSE, keeptrainfits=TRUE,
    id = NULL,     ## only used by surv.bart
    seed = 99L,    ## only used by mc.surv.bart
    mc.cores = 2L, ## ditto
    nice=19L       ## ditto
)
{
    if(.Platform$OS.type!='unix')
        stop('parallel::mcparallel/mccollect do not exist on windows')

    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    parallel::mc.reset.stream()

    if(length(y.train)==0) {
        pre <- surv.pre.bart(times, delta, x.train, x.test)

        y.train <- pre$y.train
        x.train <- pre$tx.train
        x.test  <- pre$tx.test

        if(length(binaryOffset)==0) binaryOffset <- pre$binaryOffset
    }
    else if(length(binaryOffset)==0) binaryOffset <- 0

    H <- 1
    Mx <- 2^31-1
    Nx <- max(nrow(x.train), nrow(x.test))
    
    if(Nx>Mx%/%ndpost) {
        H <- ceiling(ndpost / (Mx %/% Nx))
        ndpost <- ndpost %/% H
        ##nrow*ndpost>2Gi: due to the 2Gi limit in sendMaster
        ##(unless this limit was increased): reducing ndpost 
    }
    
    mc.cores.detected <- detectCores()

    if(mc.cores>mc.cores.detected) {
        message('The number of cores requested, ', mc.cores,
                       ',\n exceeds the number of cores detected via detectCores() ',
                       'reducing to ', mc.cores.detected)
        mc.cores <- mc.cores.detected
    }

    mc.ndpost <- ((ndpost %/% mc.cores) %/% keepevery)*keepevery

    while(mc.ndpost*mc.cores<ndpost) mc.ndpost <- mc.ndpost+keepevery

    mc.nkeep <- mc.ndpost %/% keepevery
    
    post.list <- list()

    for(h in 1:H) {
        for(i in 1:mc.cores) {
        parallel::mcparallel({psnice(value=nice);
              surv.bart(x.train=x.train, y.train=y.train, x.test=x.test,
                        k=k, power=power, base=base,
                        binaryOffset=binaryOffset,
                        ntree=ntree, numcut=numcut,
                        ndpost=mc.ndpost, nskip=nskip,
                        nkeeptrain=mc.nkeep, nkeeptest=mc.nkeep,
                        nkeeptestmean=mc.nkeep, nkeeptreedraws=mc.nkeep,
                        printevery=printevery, treesaslists=treesaslists)},
              silent=(i!=1))
              ## to avoid duplication of output
              ## capture stdout from first posterior only
        }

        post.list[[h]] <- parallel::mccollect()
    }

    if(H==1 & mc.cores==1) return(post.list[[1]][[1]])
    else {
        for(h in 1:H) for(i in mc.cores:1) {
            if(h==1 & i==mc.cores) {
                post <- post.list[[1]][[mc.cores]]
         
                p <- ncol(x.train)
                
                old.text <- paste0(as.character(mc.nkeep), ' ', as.character(ntree), ' ', as.character(p))
                old.stop <- nchar(old.text)
                
                post$treedraws$trees <- sub(old.text,
                                            paste0(as.character(H*mc.cores*mc.nkeep), ' ', as.character(ntree), ' ',
                                                   as.character(p)),
                                            post$treedraws$trees)
            }
            else {
                if(keeptrainfits) post$yhat.train <- rbind(post$yhat.train, post.list[[h]][[i]]$yhat.train)

                if(length(post$surv.train)>0)
                    post$surv.train <- rbind(post$surv.train, post.list[[h]][[i]]$surv.train)

                if(length(post$yhat.test)>0)
                    post$yhat.test <- rbind(post$yhat.test, post.list[[h]][[i]]$yhat.test)

                if(length(post$surv.test)>0)
                    post$surv.test <- rbind(post$surv.test, post.list[[h]][[i]]$surv.test)

                post$varcount <- rbind(post$varcount, post.list[[h]][[i]]$varcount)
                
                post$treedraws$trees <- paste0(post$treedraws$trees, 
                                               substr(post.list[[h]][[i]]$treedraws$trees, old.stop+2,
                                                      nchar(post.list[[h]][[i]]$treedraws$trees)))

                if(treesaslists) post$treedraws$lists <-
                                     c(post$treedraws$lists, post.list[[h]][[i]]$treedraws$lists)
                }                      

            post.list[[h]][[i]] <- NULL
            }

        if(length(post$yhat.train.mean)>0)
            post$yhat.train.mean <- apply(post$yhat.train, 2, mean)

        if(length(post$surv.train.mean)>0)
            post$surv.train.mean <- apply(post$surv.train, 2, mean)

        if(length(post$yhat.test.mean)>0)
            post$yhat.test.mean <- apply(post$yhat.test, 2, mean)

        if(length(post$surv.test.mean)>0)
            post$surv.test.mean <- apply(post$surv.test, 2, mean)

        return(post)
    }
}
