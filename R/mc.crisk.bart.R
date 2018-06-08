
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


mc.crisk.bart <- function(
    x.train = matrix(0,0,0), y.train=NULL,
    x.train2 = x.train, y.train2=NULL,
    times=NULL, delta=NULL, K=NULL,
    x.test = matrix(0,0,0), x.test2 = x.test, cond=NULL,
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
    tau.num=c(NA, 3, 6)[ntype], ##tau.num2=c(NA, 3, 6)[ntype],
    ##binaryOffset = NULL, binaryOffset2 = NULL,
    ntree = 50L, numcut = 100L,
    ndpost = 1000L, nskip = 250L,
    keepevery = 10L,
    ##nkeeptrain=ndpost, nkeeptest=ndpost,
    ##nkeeptestmean=ndpost,
    ##nkeeptreedraws=ndpost,
    printevery=100L,
    ##treesaslists=FALSE,
    ##keeptrainfits=TRUE,
    id=NULL,    ## crisk.bart only
    seed = 99L, mc.cores = 2L, nice=19L
)
{
    if(.Platform$OS.type!='unix')
        stop('parallel::mcparallel/mccollect do not exist on windows')

    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    parallel::mc.reset.stream()

    if(is.na(ntype) || ntype==1)
        stop("type argument must be set to either 'pbart' or 'lbart'")

    x.train2 <- bartModelMatrix(x.train2)
    x.test2 <- bartModelMatrix(x.test2)
    x.train <- bartModelMatrix(x.train)
    x.test <- bartModelMatrix(x.test)

    if(length(y.train)==0) {
        pre <- crisk.pre.bart(times, delta, x.train, x.test,
                              x.train2, x.test2, K=K)

        y.train <- pre$y.train
        x.train <- pre$tx.train
        x.test  <- pre$tx.test
        y.train2 <- pre$y.train2
        x.train2 <- pre$tx.train2
        x.test2  <- pre$tx.test2

        times   <- pre$times
        K       <- pre$K

        if(length(cond)==0) cond <- pre$cond
        ##if(length(binaryOffset)==0) binaryOffset <- pre$binaryOffset
        ##if(length(binaryOffset2)==0) binaryOffset2 <- pre$binaryOffset2
    }
    else {
        if(length(x.train)==0 | length(x.train2)==0)
            stop('both x.train and x.train2 must be provided')

        if(nrow(x.train)!=nrow(x.train2))
            stop('number of rows in x.train and x.train2 must be equal')

        ##if(length(binaryOffset)==0) binaryOffset <- 0
        ##if(length(binaryOffset2)==0) binaryOffset2 <- 0

        times <- unique(sort(x.train[ , 1]))
        K     <- length(times)
    }

    H <- 1
    Mx <- 2^31-1
    Nx <- 2*max(nrow(x.train), nrow(x.test))

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

    mc.ndpost <- ceiling(ndpost/mc.cores)
    ## mc.ndpost <- ((ndpost %/% mc.cores) %/% keepevery)*keepevery

    ## while(mc.ndpost*mc.cores<ndpost) mc.ndpost <- mc.ndpost+keepevery

    ## mc.nkeep <- mc.ndpost %/% keepevery

    post.list <- list()

    for(h in 1:H) {
        for(i in 1:mc.cores) {
        parallel::mcparallel({psnice(value=nice);
              crisk.bart(x.train=x.train, y.train=y.train,
                         x.train2=x.train2, y.train2=y.train2,
                         x.test=x.test, x.test2=x.test2, cond=cond,
                         sparse=sparse, theta=theta, omega=omega,
                         a=a, b=b, augment=augment, rho=rho, rho2=rho2,
                         xinfo=xinfo, xinfo2=xinfo2, usequants=usequants,
                         ##cont=cont,
                         rm.const=rm.const, type=type,
                         k=k, power=power, base=base,
                         ##binaryOffset=binaryOffset,
                         ##binaryOffset2=binaryOffset2,
                         offset=offset, offset2=offset2,
                         tau.num=tau.num, ##tau.num2=tau.num2,
                         ntree=ntree, numcut=numcut,
                         ndpost=mc.ndpost, nskip=nskip,
                         keepevery = keepevery,
                         ##nkeeptrain=mc.ndpost, nkeeptest=mc.ndpost,
                         ##nkeeptestmean=mc.ndpost,
                         ##nkeeptreedraws=mc.ndpost,
                         printevery=printevery)},
                         ##treesaslists=FALSE,
                         ##keeptrainfits=TRUE)},
              silent=(i!=1))
              ## to avoid duplication of output
              ## capture stdout from first posterior only
        }

        post.list[[h]] <- parallel::mccollect()
    }

    if((H==1 & mc.cores==1) | attr(post.list[[1]][[1]], 'class')!='criskbart')
        return(post.list[[1]][[1]])
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

                p <- ncol(x.train2[ , post$rm.const2])
                old.text <- paste0(as.character(mc.ndpost), ' ',
                                   as.character(ntree), ' ', as.character(p))
                old.stop2 <- nchar(old.text)

                post$treedraws2$trees <- sub(old.text,
                                            paste0(as.character(post$ndpost),
                                                   ' ', as.character(ntree),
                                                   ' ', as.character(p)),
                                            post$treedraws2$trees)
            }
            else {
                if(length(x.test)>0) {
                    post$yhat.test <- rbind(post$yhat.test,
                                            post.list[[h]][[i]]$yhat.test)
                    post$yhat.test2 <- rbind(post$yhat.test2,
                                             post.list[[h]][[i]]$yhat.test2)
                    post$prob.test <- rbind(post$prob.test,
                                            post.list[[h]][[i]]$prob.test)
                    post$prob.test2 <- rbind(post$prob.test2,
                                             post.list[[h]][[i]]$prob.test2)
                    post$cif.test <- rbind(post$cif.test,
                                           post.list[[h]][[i]]$cif.test)
                    post$cif.test2 <- rbind(post$cif.test2,
                                            post.list[[h]][[i]]$cif.test2)
                    post$surv.test <- rbind(post$surv.test,
                                            post.list[[h]][[i]]$surv.test)
                }

                post$varcount <- rbind(post$varcount,
                                       post.list[[h]][[i]]$varcount)
                post$varcount2 <- rbind(post$varcount2,
                                        post.list[[h]][[i]]$varcount2)
                post$varprob <- rbind(post$varprob,
                                      post.list[[h]][[i]]$varprob)
                post$varprob2 <- rbind(post$varprob2,
                                       post.list[[h]][[i]]$varprob2)

                post$treedraws$trees <- paste0(post$treedraws$trees,
                                               substr(post.list[[h]][[i]]$treedraws$trees, old.stop+2,
                                                      nchar(post.list[[h]][[i]]$treedraws$trees)))

                post$treedraws2$trees <- paste0(post$treedraws2$trees,
                                               substr(post.list[[h]][[i]]$treedraws2$trees, old.stop2+2,
                                                      nchar(post.list[[h]][[i]]$treedraws2$trees)))

                ## if(treesaslists) {
                ##     post$treedraws$lists <-
                ##         c(post$treedraws$lists, post.list[[h]][[i]]$treedraws$lists)

                ##     post$treedraws2$lists <-
                ##         c(post$treedraws2$lists, post.list[[h]][[i]]$treedraws2$lists)
                ## }
            }

            post.list[[h]][[i]] <- NULL
            }

        if(length(x.test)>0) {
            post$prob.test.mean <- apply(post$prob.test, 2, mean)
            post$prob.test2.mean <- apply(post$prob.test2, 2, mean)
            post$cif.test.mean <- apply(post$cif.test, 2, mean)
            post$cif.test2.mean <- apply(post$cif.test2, 2, mean)
            post$surv.test.mean <- apply(post$surv.test, 2, mean)
        }

        post$varcount.mean <- apply(post$varcount, 2, mean)
        post$varcount2.mean <- apply(post$varcount2, 2, mean)
        post$varprob.mean <- apply(post$varprob, 2, mean)
        post$varprob2.mean <- apply(post$varprob2, 2, mean)

        attr(post, 'class') <- 'criskbart'

        return(post)
    }
}
