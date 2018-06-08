
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

mc.recur.bart <- function(
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
    if(.Platform$OS.type!='unix')
        stop('parallel::mcparallel/mccollect do not exist on windows')

    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    parallel::mc.reset.stream()

    if(is.na(ntype) || ntype==1)
        stop("type argument must be set to either 'pbart' or 'lbart'")

    x.train <- bartModelMatrix(x.train)
    ##x.test <- bartModelMatrix(x.test)

    if(length(y.train)==0) {
        recur <- recur.pre.bart(times, delta, x.train)
        ##recur <- recur.pre.bart(times, delta, x.train, x.test)

        y.train <- recur$y.train
        x.train <- recur$tx.train
        x.test  <- recur$tx.test
    }
    ##     if(length(binaryOffset)==0) {
    ##         lambda <- sum(delta, na.rm=TRUE)/
    ##             sum(apply(times, 1, max, na.rm=TRUE))
    ##         ##lambda <- sum(delta)/sum(times[ , ncol(times)])
    ##         binaryOffset <- qnorm(1-exp(-lambda))
    ##     }
    ## }
    ## else if(length(binaryOffset)==0) binaryOffset <- 0

    H <- 1
    Mx <- 2^31-1
    Nx <- max(nrow(x.train), nrow(x.test))

    if(Nx>Mx%/%ndpost) {
        H <- ceiling(ndpost / (Mx %/% Nx))
        ndpost <- ndpost %/% H
        ##nrow*ndpost>2Gi!
        ##due to the 2Gi limit in sendMaster, breaking run into H parts
        ##this bug/feature might be addressed in R-devel post 3.3.2
        ##i.e., the fix is NOT in R version 3.3.2
        ##will revisit once there is a fix in an official R release
        ##New Features entry for R-devel post 3.3.2
## The unexported low-level functions in package parallel for passing
## serialized R objects to and from forked children now support long
## vectors on 64-bit platforms. This removes some limits on
## higher-level functions such as mclapply()
    }

    mc.cores.detected <- detectCores()

    if(mc.cores>mc.cores.detected) mc.cores->mc.cores.detected
        ## warning(paste0('The number of cores requested, mc.cores=', mc.cores,
        ##                ',\n exceeds the number of cores detected via detectCores() ',
        ##                'which yields ', mc.cores.detected, ' .'))

    mc.ndpost <- ceiling(ndpost/mc.cores)
    ## mc.ndpost <- ((ndpost %/% mc.cores) %/% keepevery)*keepevery

    ## while(mc.ndpost*mc.cores<ndpost) mc.ndpost <- mc.ndpost+keepevery

    ## mc.nkeep <- mc.ndpost %/% keepevery

    post.list <- list()

    for(h in 1:H) {
        for(i in 1:mc.cores) {
            parallel::mcparallel({psnice(value=nice);
                recur.bart(x.train=x.train, y.train=y.train,
                           x.test=x.test, x.test.nogrid=x.test.nogrid,
                           sparse=sparse, theta=theta, omega=omega,
                           a=a, b=b, augment=augment, rho=rho,
                           xinfo=xinfo, usequants=usequants,
                           ##cont=cont,
                           rm.const=rm.const, type=type,
                           k=k, power=power, base=base,
                           offset=offset, tau.num=tau.num,
                           ##binaryOffset=binaryOffset,
                           ntree=ntree, numcut=numcut,
                           ndpost=mc.ndpost, nskip=nskip, keepevery=keepevery,
                           ##nkeeptrain=mc.ndpost, nkeeptest=mc.ndpost,
                           ##nkeeptestmean=mc.ndpost,
                           ##nkeeptreedraws=mc.ndpost,
                           printevery=printevery, ##treesaslists=treesaslists,
                           keeptrainfits=keeptrainfits)},
                silent=(i!=1))
            ## to avoid duplication of output
            ## capture stdout from first posterior only
        }

        post.list[[h]] <- parallel::mccollect()
    }

    if((H==1 & mc.cores==1) | attr(post.list[[1]][[1]], 'class')!='recurbart') return(post.list[[1]][[1]])
    else {
        for(h in 1:H) for(i in mc.cores:1) {
            if(h==1 & i==mc.cores) {
                post <- post.list[[1]][[mc.cores]]
                post$ndpost <- H*mc.cores*mc.ndpost
                p <- ncol(x.train[ , post$rm.const])

                old.text <- paste0(as.character(mc.ndpost), ' ',
                                   as.character(ntree), ' ', as.character(p))
                old.stop <- nchar(old.text)

                post$treedraws$trees <- sub(old.text,
                                            paste0(as.character(post$ndpost),
                                                   ' ', as.character(ntree),
                                                   ' ', as.character(p)),
                                            post$treedraws$trees)
            }
            else {
                if(keeptrainfits) {
                    post$yhat.train <- rbind(post$yhat.train,
                                             post.list[[h]][[i]]$yhat.train)
                    post$prob.train <- rbind(post$prob.train,
                                             post.list[[h]][[i]]$prob.train)
                    post$haz.train <- rbind(post$haz.train,
                                            post.list[[h]][[i]]$haz.train)
                    post$cum.train <- rbind(post$cum.train,
                                            post.list[[h]][[i]]$cum.train)
                }

                if(length(x.test)>0) {
                    post$yhat.test <- rbind(post$yhat.test, post.list[[h]][[i]]$yhat.test)
                    post$prob.test <- rbind(post$prob.test, post.list[[h]][[i]]$prob.test)
                    post$haz.test <- rbind(post$haz.test, post.list[[h]][[i]]$haz.test)
                    if(!x.test.nogrid) post$cum.test <- rbind(post$cum.test, post.list[[h]][[i]]$cum.test)
                }

                ## if(length(post$sigma)>0)
                ##     post$sigma <- c(post$sigma, post.list[[h]][[i]]$sigma)

                post$varcount <- rbind(post$varcount, post.list[[h]][[i]]$varcount)
                post$varprob <- rbind(post$varprob, post.list[[h]][[i]]$varprob)

                post$treedraws$trees <- paste0(post$treedraws$trees,
                                               substr(post.list[[h]][[i]]$treedraws$trees, old.stop+2,
                                                      nchar(post.list[[h]][[i]]$treedraws$trees)))

                ## if(treesaslists) post$treedraws$lists <-
                ##                      c(post$treedraws$lists, post.list[[h]][[i]]$treedraws$lists)

            }

            post.list[[h]][[i]] <- NULL
        }

        ##post$yhat.train.mean <- apply(post$yhat.train, 2, mean)

        if(keeptrainfits) {
            post$prob.train.mean <- apply(post$prob.train, 2, mean)
            post$haz.train.mean <- apply(post$haz.train, 2, mean)
            post$cum.train.mean <- apply(post$cum.train, 2, mean)
        }

        if(length(x.test)>0) {
            post$prob.test.mean <- apply(post$prob.test, 2, mean)
            post$haz.test.mean <- apply(post$haz.test, 2, mean)
            if(!x.test.nogrid) post$cum.test.mean <- apply(post$cum.test, 2, mean)
        }

        post$varcount.mean <- apply(post$varcount, 2, mean)
        post$varprob.mean <- apply(post$varprob, 2, mean)

        attr(post, 'class') <- 'recurbart'

        return(post)
    }
}
