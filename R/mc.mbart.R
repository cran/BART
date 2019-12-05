
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


mc.mbart <- function(
                     x.train, y.train,
                     x.test=matrix(0,0,0), type='pbart',
                     ntype=as.integer(
                         factor(type,
                                levels=c('wbart', 'pbart', 'lbart'))),
                     sparse=FALSE, theta=0, omega=1,
                     a=0.5, b=1, augment=FALSE, rho=NULL,
                     xinfo=matrix(0,0,0), usequants=FALSE,
                     rm.const=TRUE,
                     k=2, power=2, base=0.95,
                     ##sigest=NA, sigdf=3, sigquant=0.90, lambda=NA,
                     tau.num=c(NA, 3, 6)[ntype],
                     offset=NULL, ##w=rep(1, length(y.train)),
                     ntree=c(200L, 50L, 50L)[ntype], numcut=100L,
                     ndpost=1000L, nskip=100L,
                     keepevery=c(1L, 10L, 10L)[ntype],
                     printevery=100L, transposed=FALSE,
                     hostname=FALSE,
                     mc.cores = 2L, nice = 19L, seed = 99L
                     )
{
    if(.Platform$OS.type!='unix')
        stop('parallel::mcparallel/mccollect do not exist on windows')

    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    parallel::mc.reset.stream()

    if(type=='wbart' || is.na(ntype))
        stop("type argument must be set to either 'pbart' or 'lbart'")

    cats <- unique(sort(y.train))
    K <- length(cats)
    if(K<2)
        stop("there must be at least 2 categories")

    L <- length(offset)
    if(!(L %in% c(0, K)))
        stop(paste0("length of offset argument must be 0 or ", K))

    L <- K-1
    
    if(!transposed) {
        temp = bartModelMatrix(x.train, numcut, usequants=usequants,
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
    }

    mc.cores.detected <- detectCores()

    if(mc.cores>mc.cores.detected) mc.cores <- mc.cores.detected

    mc.ndpost <- ceiling(ndpost/mc.cores)

    for(i in 1:mc.cores) {
        parallel::mcparallel({psnice(value=nice);
            mbart(x.train=x.train, y.train=y.train,
                  x.test=x.test,
                  type=type, ntype=ntype,
                  sparse=sparse, theta=theta, omega=omega,
                  a=a, b=b, augment=augment, rho=rho,
                  xinfo=xinfo, usequants=usequants,
                  rm.const=rm.const,
                  k=k, power=power, base=base,
                  tau.num=tau.num,
                  offset=offset,
                  ntree=ntree, numcut=numcut,
                  ndpost=mc.ndpost, nskip=nskip,
                  keepevery=keepevery,
                  printevery=printevery, transposed=TRUE,
                  hostname=hostname)},
            silent=(i!=1))
        ## to avoid duplication of output
        ## capture stdout from first posterior only
    }

    post.list <- parallel::mccollect()

    post <- post.list[[1]]

    if(mc.cores==1 | attr(post, 'class')!='mbart') return(post)
    else {
        if(class(rm.const)[1]!='logical') post$rm.const <- rm.const
        post$ndpost <- mc.cores*mc.ndpost
        p <- nrow(x.train[ , post$rm.const])

        if(length(rm.const)==0) rm.const <- 1:p

        post$rm.const <- rm.const

        old.text <- paste0(as.character(mc.ndpost), ' ', as.character(ntree),
                           ' ', as.character(p))
        old.stop <- nchar(old.text)

        for(j in 1:L)
            post$treedraws$trees[[j]] <- sub(old.text,
                                       paste0(as.character(post$ndpost), ' ',
                                              as.character(ntree), ' ',
                                              as.character(p)),
                                       post$treedraws$trees[[j]])

        keeptestfits <- length(x.test)>0

        for(i in 2:mc.cores) {

            if(keeptestfits) {
                post$yhat.test <- rbind(post$yhat.test,
                                        post.list[[i]]$yhat.test)
                post$prob.test <- rbind(post$prob.test,
                                        post.list[[i]]$prob.test)
            }

            for(j in 1:L) {
                post$varcount[[j]] <- rbind(post$varcount[[j]],
                                            post.list[[i]]$varcount[[j]])
                post$varprob[[j]] <- rbind(post$varprob[[j]],
                                           post.list[[i]]$varprob[[j]])

                post$treedraws$trees[[j]] <- paste0(post$treedraws$trees[[j]],
                         substr(post.list[[i]]$treedraws$trees[[j]], old.stop+2,
                                    nchar(post.list[[i]]$treedraws$trees[[j]])))
            }

        }

        ##post$prob.train.mean <- apply(post$prob.train, 2, mean)

        if(keeptestfits) post$prob.test.mean <- apply(post$prob.test, 2, mean)

        for(j in 1:L) {
            post$varcount.mean[j, ] <- apply(post$varcount[[j]], 2, mean)
            post$varprob.mean[j, ] <- apply(post$varprob[[j]], 2, mean)
        }
        
        attr(post, 'class') <- 'mbart'

        return(post)
    }
}
