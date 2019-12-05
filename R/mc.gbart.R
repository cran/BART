
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


mc.gbart <- function(
                     x.train, y.train,
                     x.test=matrix(0,0,0), type='wbart',
                     ntype=as.integer(
                         factor(type,
                                levels=c('wbart', 'pbart', 'lbart'))),
                     sparse=FALSE, theta=0, omega=1,
                     a=0.5, b=1, augment=FALSE, rho=NULL,
                     xinfo=matrix(0,0,0), usequants=FALSE,
                     rm.const=TRUE,
                     sigest=NA, sigdf=3, sigquant=0.90,
                     k=2, power=2, base=0.95,
                     ##sigmaf=NA,
                     lambda=NA, tau.num=c(NA, 3, 6)[ntype],
                     ##tau.interval=0.9973,
                     offset=NULL, w=rep(1, length(y.train)),
                     ntree=c(200L, 50L, 50L)[ntype], numcut=100L,
                     ndpost=1000L, nskip=100L,
                     keepevery=c(1L, 10L, 10L)[ntype],
                     printevery=100L, transposed=FALSE,
                     hostname=FALSE,
                     mc.cores = 2L, nice = 19L, seed = 99L
                     )
{
    if(is.na(ntype))
        stop("type argument must be set to either 'wbart', 'pbart' or 'lbart'")

    check <- unique(sort(y.train))

    if(length(check)==2) {
        if(!all(check==0:1))
            stop('Binary y.train must be coded as 0 and 1')
        if(type=='wbart')
            stop("The outcome is binary so set type to 'pbart' or 'lbart'")
    }

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
            gbart(x.train=x.train, y.train=y.train,
                  x.test=x.test, type=type, ntype=ntype,
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
                  transposed=TRUE, hostname=hostname)},
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

        if(type=='wbart') sigma <- post$sigma[-(1:nskip)]
        
        for(i in 2:mc.cores) {
            post$hostname[i] <- post.list[[i]]$hostname
            
            post$yhat.train <- rbind(post$yhat.train,
                                     post.list[[i]]$yhat.train)

            if(keeptest) post$yhat.test <- rbind(post$yhat.test,
                                                 post.list[[i]]$yhat.test)

            if(type=='wbart') {
                post$sigma <- cbind(post$sigma, post.list[[i]]$sigma)
                sigma <- c(sigma, post.list[[i]]$sigma[-(1:nskip)])
            }

            post$varcount <- rbind(post$varcount, post.list[[i]]$varcount)
            post$varprob <- rbind(post$varprob, post.list[[i]]$varprob)

            post$treedraws$trees <-
                paste0(post$treedraws$trees,
                       substr(post.list[[i]]$treedraws$trees, old.stop+2,
                              nchar(post.list[[i]]$treedraws$trees)))

            post$proc.time['elapsed'] <-
                max(post$proc.time['elapsed'],
                    post.list[[i]]$proc.time['elapsed'])

            for(j in 1:5)
                if(j!=3)
                    post$proc.time[j] <-
                        post$proc.time[j]+post.list[[i]]$proc.time[j]
        }
        
        n=length(y.train)
        Y=t(matrix(y.train, nrow=n, ncol=post$ndpost))
        
        if(type=='wbart') {
            post$yhat.train.mean <- apply(post$yhat.train, 2, mean)
            SD=matrix(sigma, nrow=post$ndpost, ncol=n)
            ##CPO=1/apply(1/dnorm(Y, post$yhat.train, SD), 2, mean)
            log.pdf=dnorm(Y, post$yhat.train, SD, TRUE)
            post$sigma.mean=mean(SD[ , 1])
            
            if(keeptest)
                post$yhat.test.mean <- apply(post$yhat.test, 2, mean)
        } else {
            post$prob.train.mean <- apply(post$prob.train, 2, mean)
            ##CPO=1/apply(1/dbinom(Y, 1, post$prob.train), 2, mean)
            log.pdf=dbinom(Y, 1, post$prob.train, TRUE)

            if(keeptest)
                post$prob.test.mean <- apply(post$prob.test, 2, mean)
        }

        min.log.pdf=t(matrix(apply(log.pdf, 2, min),
                             nrow=n, ncol=post$ndpost))
        log.CPO=log(post$ndpost)+min.log.pdf[1, ]-
            log(apply(exp(min.log.pdf-log.pdf), 2, sum))

        post$LPML=sum(log.CPO) 
        ##post$LPML=sum(log(CPO))

        post$varcount.mean <- apply(post$varcount, 2, mean)
        post$varprob.mean <- apply(post$varprob, 2, mean)

        attr(post, 'class') <- type

        return(post)
    }
}
