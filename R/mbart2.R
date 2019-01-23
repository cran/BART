
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


mbart2 <- function(
                   x.train, y.train,
                   x.test=matrix(0,0,0), type='lbart',
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
    if(type=='wbart' || is.na(ntype))
        stop("type argument must be set to either 'pbart' or 'lbart'")

    cats <- unique(sort(y.train))
    K <- length(cats)

    if(K<2)
        stop("there must be at least 2 categories")

    L=length(offset)

    if(!(L %in% c(0, K)))
        stop(paste0("length of offset argument must be 0 or ", K))

    if(!transposed) {
        temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                               xinfo=xinfo, rm.const=rm.const)
        x.train = t(temp$X)
        numcut = temp$numcut
        xinfo = temp$xinfo
        if(length(x.test)>0) {
            x.test = bartModelMatrix(x.test)
            x.test = t(x.test[ , temp$rm.const])
        }
        rm.const <- temp$rm.const
        rm(temp)
    }

    post <- list()
    post$K <- K
    post$cats <- cats
    N <- length(y.train)
    P <- nrow(x.train) ## transposed
    post$yhat.train <- matrix(nrow=ndpost, ncol=N*K)
    post$prob.train <- matrix(nrow=ndpost, ncol=N*K)
    post$tot.train <- matrix(0, nrow=ndpost, ncol=N)
    if(length(x.test)) {
        Q <- ncol(x.test) ## transposed
        post$yhat.test <- matrix(nrow=ndpost, ncol=Q*K)
        post$prob.test <- matrix(nrow=ndpost, ncol=Q*K)
        post$tot.test <- matrix(0, nrow=ndpost, ncol=Q)
    }
    else Q <- 0

    post$varcount <- as.list(1:K)
    post$varprob <- as.list(1:K)
    post$varcount.mean <- matrix(nrow=K, ncol=P)
    post$varprob.mean <- matrix(nrow=K, ncol=P)
    post$offset <- 0
    post$treedraws <- list()
    post$treedraws$trees <- as.list(1:K)
    post.list <- as.list(1:K)

    for(h in 1:K) {
        post.list[[h]] <-
            gbart(x.train=x.train,
                  y.train=(y.train==h)*1,
                  x.test=x.test,
                  type=type, ntype=ntype,
                  sparse=sparse, theta=theta, omega=omega,
                  a=a, b=b, augment=augment, rho=rho,
                  xinfo=xinfo, usequants=usequants,
                  rm.const=rm.const,
                  k=k, power=power, base=base,
                  tau.num=tau.num,
                  offset=offset[h],
                  ntree=ntree, numcut=numcut,
                  ndpost=ndpost, nskip=nskip,
                  keepevery=keepevery,
                  printevery=printevery, transposed=TRUE,
                  hostname=hostname)

        if(attr(post.list[[h]], 'class')!=type) return(post.list[[h]])

        for(i in 1:N) {
            j <- (i-1)*K+h
            post$yhat.train[ , j] <- post.list[[h]]$yhat.train[ , i]
            post$tot.train[ , i] <-
                post$tot.train[ , i]+exp(post$yhat.train[ , j])
            if(h==K) 
                for(j in (i-1)*K+1:K)
                    post$prob.train[ , j] <-
                        exp(post$yhat.train[ , j])/post$tot.train[ , i]
        }

        if(Q>0) {
            for(i in 1:Q) {
                j <- (i-1)*K+h
                post$yhat.test[ , j] <- post.list[[h]]$yhat.test[ , i]
                post$tot.test[ , i] <-
                    post$tot.test[ , i]+exp(post$yhat.test[ , j])
                if(h==K) 
                    for(j in (i-1)*K+1:K)
                        post$prob.test[ , j] <-
                            exp(post$yhat.test[ , j])/post$tot.test[ , i]
            }
        }

        post$varcount[[h]] <- post.list[[h]]$varcount
        post$varprob[[h]] <- post.list[[h]]$varprob

        for(j in 1:P) {
            ##post$varcount[ , j, h] <- post.list[[h]]$varcount[ , j]
            post$varcount.mean[h, j] <- post.list[[h]]$varcount.mean[j]
            ##post$varprob[ , j, h] <- post.list[[h]]$varprob[ , j]
            post$varprob.mean[h, j] <- post.list[[h]]$varprob.mean[j]
        }

        post$offset[h] <- post.list[[h]]$offset
        ##post$rm.const[[h]] <- post.list[[h]]$rm.const
        post$treedraws$trees[[h]] <- post.list[[h]]$treedraws$trees
    }

    post$treedraws$cutpoints <- post.list[[1]]$treedraws$cutpoints
    dimnames(post$varcount.mean)[[2]] <- dimnames(post$varcount[[1]])[[2]]
    dimnames(post$varprob.mean)[[2]] <- dimnames(post$varprob[[1]])[[2]]
    post$rm.const <- post.list[[1]]$rm.const
    post$type <- type

    post$prob.train.mean <- apply(post$prob.train, 2, mean)
    if(Q>0) {
        post$prob.test.mean <- apply(post$prob.test, 2, mean)
        ##post$comp.test.mean <- apply(post$comp.test, 2, mean)
    }

    post$tot.train <- NULL
    post$tot.test <- NULL
    attr(post, 'class') <- 'mbart2'

    return(post)
}
