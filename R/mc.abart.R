
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


mc.abart <- function(
                     x.train, times, delta,
                     x.test=matrix(0,0,0), K=100,
                     type='abart', ntype=1,
                     sparse=FALSE, theta=0, omega=1,
                     a=0.5, b=1, augment=FALSE, rho=NULL,
                     xinfo=matrix(0,0,0), usequants=FALSE,
                     rm.const=TRUE,
                     sigest=NA, sigdf=3, sigquant=0.90,
                     k=2, power=2, base=0.95,
                     ##sigmaf=NA,
                     lambda=NA, tau.num=c(NA, 3, 6)[ntype],
                     ##tau.interval=0.9973,
                     offset=NULL, w=rep(1, length(times)),
                     ntree=c(200L, 50L, 50L)[ntype], numcut=100L,
                     ndpost=1000L, nskip=100L,
                     keepevery=c(1L, 10L, 10L)[ntype],
                     printevery=100L, transposed=FALSE,
                     mc.cores = 2L, nice = 19L, seed = 99L
                     )
{
    if(type!='abart') stop('type must be "abart"')
    if(ntype!=1) stop('ntype must be 1')

    if(.Platform$OS.type!='unix')
        stop('parallel::mcparallel/mccollect do not exist on windows')

    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    parallel::mc.reset.stream()

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
            abart(x.train=x.train, times=times, delta=delta,
                  x.test=x.test, K=K, type=type, ntype=ntype,
                  sparse=sparse, theta=theta, omega=omega,
                  a=a, b=b, augment=augment, rho=rho,
                  xinfo=xinfo, usequants=usequants,
                  rm.const=rm.const,
                  sigest=sigest, sigdf=sigdf, sigquant=sigquant,
                  k=k, power=power, base=base,
                  ##sigmaf=sigmaf,
                  lambda=lambda, tau.num=tau.num,
                  ##tau.interval=tau.interval,
                  offset=offset,
                  w=w, ntree=ntree, numcut=numcut,
                  ndpost=mc.ndpost, nskip=nskip,
                  keepevery=keepevery, printevery=printevery,
                  transposed=TRUE)},
            silent=(i!=1))
        ## to avoid duplication of output
        ## capture stdout from first posterior only
    }

    post.list <- parallel::mccollect()

    post <- post.list[[1]]

    if(mc.cores==1 | attr(post, 'class')!=type) return(post)
    else {
        if(class(rm.const)[1]!='logical') post$rm.const <- rm.const

        post$ndpost <- mc.cores*mc.ndpost

        p <- nrow(x.train[post$rm.const, ])

        old.text <- paste0(as.character(mc.ndpost), ' ', as.character(ntree),
                           ' ', as.character(p))
        old.stop <- nchar(old.text)

        post$treedraws$trees <- sub(old.text,
                                    paste0(as.character(post$ndpost), ' ',
                                           as.character(ntree), ' ',
                                           as.character(p)),
                                    post$treedraws$trees)

        keeptest <- length(x.test)>0

        for(i in 2:mc.cores) {
            post$sigma <- cbind(post$sigma, post.list[[i]]$sigma)

            post$yhat.train <- rbind(post$yhat.train,
                                     post.list[[i]]$yhat.train)

            post$surv.train <- rbind(post$surv.train,
                                     post.list[[i]]$surv.train)

            if(keeptest) {
                post$yhat.test <- rbind(post$yhat.test,
                                        post.list[[i]]$yhat.test)

                post$surv.test <- rbind(post$surv.test,
                                        post.list[[i]]$surv.test)
            }

            post$varcount <- rbind(post$varcount, post.list[[i]]$varcount)
            post$varprob <- rbind(post$varprob, post.list[[i]]$varprob)

            post$treedraws$trees <- paste0(post$treedraws$trees,
                                           substr(post.list[[i]]$treedraws$trees, old.stop+2,
                                                  nchar(post.list[[i]]$treedraws$trees)))

            post$proc.time['elapsed'] <- max(post$proc.time['elapsed'],
                                             post.list[[i]]$proc.time['elapsed'])
            for(j in 1:5)
                if(j!=3)
                    post$proc.time[j] <- post$proc.time[j]+post.list[[i]]$proc.time[j]
        }

        if(type=='abart') {
            post$yhat.train.mean <- apply(post$yhat.train, 2, mean)
            post$surv.train.mean <- apply(post$surv.train, 2, mean)

            if(keeptest) {
                post$yhat.test.mean <- apply(post$yhat.test, 2, mean)
                post$surv.test.mean <- apply(post$surv.test, 2, mean)
            }
        } else {
            post$prob.train.mean <- apply(post$prob.train, 2, mean)

            if(keeptest)
                post$prob.test.mean <- apply(post$prob.test, 2, mean)
        }

        post$varcount.mean <- apply(post$varcount, 2, mean)
        post$varprob.mean <- apply(post$varprob, 2, mean)

        attr(post, 'class') <- type

        return(post)
    }
}
