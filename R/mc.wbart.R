
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

mc.wbart <- function(
    x.train, y.train, x.test=matrix(0.0,0,0),
    sparse=FALSE, theta=0, omega=1,
    a=0.5, b=1, augment=FALSE, rho=NULL,
    xinfo=matrix(0.0,0,0), usequants=FALSE,
    cont=FALSE, rm.const=TRUE,
    sigest=NA, sigdf=3, sigquant=0.90,
    k=2.0, power=2.0, base=.95,
    sigmaf=NA, lambda=NA,
    fmean=mean(y.train),
    w=rep(1,length(y.train)),
    ntree=200L, numcut=100L,
    ndpost=1000L, nskip=100L,
    keepevery=1L, printevery=100L,
    keeptrainfits=TRUE, transposed=FALSE,
    ##treesaslists=FALSE,
    mc.cores = 2L, nice = 19L,
    seed = 99L
    )
{
    if(.Platform$OS.type!='unix')
        stop('parallel::mcparallel/mccollect do not exist on windows')

    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    parallel::mc.reset.stream()

    if(!transposed) {
        temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                               cont=cont, xinfo=xinfo, rm.const=rm.const)
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
        ## warning(paste0('The number of cores requested, mc.cores=', mc.cores,
        ##                ',\n exceeds the number of cores detected via detectCores() ',
        ##                'which yields ', mc.cores.detected, ' .'))

    mc.ndpost <- ceiling(ndpost/mc.cores)

    for(i in 1:mc.cores) {
        parallel::mcparallel({psnice(value=nice);
                   wbart(x.train=x.train, y.train=y.train, x.test=x.test,
                         sparse=sparse, theta=theta, omega=omega,
                         a=a, b=b, augment=augment, rho=rho,
                         xinfo=xinfo,
                         sigest=sigest, sigdf=sigdf, sigquant=sigquant,
                         k=k, power=power, base=base,
                         sigmaf=sigmaf, lambda=lambda, fmean=fmean, w=w,
                         ntree=ntree, numcut=numcut,
                         ndpost=mc.ndpost, nskip=nskip, keepevery=keepevery,
                         printevery=printevery, transposed=TRUE)},
                         ##treesaslists=treesaslists)},
                   silent=(i!=1))
                   ## to avoid duplication of output
                   ## capture stdout from first posterior only
    }

    post.list <- parallel::mccollect()

    post <- post.list[[1]]

    ## sigma.len <- length(post$sigma)
    ## if(sigma.len>mc.ndpost) {
    ##     sigma.beg <- 1+sigma.len-mc.ndpost
    ##     post$sigma <- post$sigma[sigma.beg:sigma.len]
    ## }

    if(mc.cores==1 | attr(post, 'class')!='wbart') return(post)
    else {
        if(class(rm.const)[1]!='logical') post$rm.const <- rm.const
        post$ndpost <- mc.cores*mc.ndpost
        p <- nrow(x.train[post$rm.const, ])
        ##p <- nrow(x.train[ , post$rm.const])

        ## if(length(rm.const)==0) rm.const <- 1:p
        ## post$rm.const <- rm.const

        old.text <- paste0(as.character(mc.ndpost), ' ', as.character(ntree),
                           ' ', as.character(p))
        old.stop <- nchar(old.text)

        post$treedraws$trees <- sub(old.text,
                                    paste0(as.character(post$ndpost), ' ',
                                           as.character(ntree), ' ',
                                           as.character(p)),
                                    post$treedraws$trees)

        for(i in 2:mc.cores) {
            post$yhat.train <- rbind(post$yhat.train,
                                     post.list[[i]]$yhat.train)

            if(length(post$yhat.test)>0)
                post$yhat.test <- rbind(post$yhat.test,
                                        post.list[[i]]$yhat.test)

            ## if(sigma.len>0)
            ##     post$sigma <- c(post$sigma, post.list[[i]]$sigma[sigma.beg:sigma.len])

            post$sigma <- cbind(post$sigma, post.list[[i]]$sigma)

            post$treedraws$trees <- paste0(post$treedraws$trees,
                                           substr(post.list[[i]]$treedraws$trees, old.stop+2,
                                                  nchar(post.list[[i]]$treedraws$trees)))

            ## if(treesaslists) post$treedraws$lists <-
            ##                      c(post$treedraws$lists, post.list[[i]]$treedraws$lists)

            if(length(post$varcount)>0) {
                post$varcount <- rbind(post$varcount, post.list[[i]]$varcount)
                post$varprob <- rbind(post$varprob, post.list[[i]]$varprob)
            }

            post$proc.time['elapsed'] <- max(post$proc.time['elapsed'],
                                                 post.list[[i]]$proc.time['elapsed'])
            for(j in 1:5)
                if(j!=3)
                    post$proc.time[j] <- post$proc.time[j]+post.list[[i]]$proc.time[j]
        }

        if(length(post$yhat.train.mean)>0)
            post$yhat.train.mean <- apply(post$yhat.train, 2, mean)

        if(length(post$yhat.test.mean)>0)
            post$yhat.test.mean <- apply(post$yhat.test, 2, mean)

        if(length(post$varcount)>0) {
            post$varcount.mean <- apply(post$varcount, 2, mean)
            post$varprob.mean <- apply(post$varprob, 2, mean)
        }

        attr(post, 'class') <- 'wbart'

        return(post)
    }
}
