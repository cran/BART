
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


rs.wbart <- function(
    x.train, y.train, x.test = matrix(0.0, 0L, 0L),
    C = floor(length(y.train)/2000),
    k = 2.0, ## BEWARE: do NOT use k for other purposes below
    power = 2.0, base = 0.95,
    ntree=50L, numcut=100L,
    ndpost=1000L, nskip=100L,
    keepevery=1L, printevery=100L,
    keeptrainfits=FALSE,
    transposed=FALSE,
    mc.cores = 2L, nice = 19L,
    seed = 99L
)
{
    if(C<=1) stop('The number of shards must be >1')

    if(.Platform$OS.type!='unix')
        stop('parallel::mcparallel/mccollect do not exist on windows')

    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    parallel::mc.reset.stream()

    if(!transposed) {
        temp = bartModelMatrix(x.train, numcut)
        x.train = t(temp$X)
        numcut = temp$numcut
        xinfo = temp$xinfo
        rm(temp)
    }

    mc.cores.detected <- detectCores()

    if(mc.cores>mc.cores.detected) mc.cores <- mc.cores.detected

    mc.ndpost <- ceiling(ndpost/mc.cores)

    p <- nrow(x.train)

    rs <- stratrs(y.train, C)

    post.list <- as.list(1:C)
    shard <- as.list(1:C)

    for(h in 1:C) {
        shard[[h]] <- rs==h
        for(i in 1:mc.cores) {
            parallel::mcparallel({psnice(value=nice);
                wbart(x.train=x.train[ , shard[[h]] ],
                      y.train=y.train[shard[[h]] ],
                      xinfo=xinfo,
                      k=k, power=power, base=base,
                      ntree=ntree, numcut=numcut,
                      ndpost=mc.ndpost, nskip=nskip, keepevery=keepevery,
                      printevery=printevery, transposed=TRUE)},
                silent=(i!=1))
            ## to avoid duplication of output
            ## capture stdout from first posterior only
        }

        post.list[[h]] <- parallel::mccollect()

        if(h==1) {
            post <- post.list[[1]][[1]]

            if(attr(post, 'class')!='wbart') return(post)

            post$yhat.shard <- matrix(nrow=mc.ndpost*mc.cores,
                                      ncol=length(y.train))
            post$yhat.shard[1:mc.ndpost, shard[[1]] ] <- post$yhat.train
            post$yhat.train <- NULL
            post$yhat.train.mean <- NULL
        }
    }

    old.text <- paste0(as.character(mc.ndpost), ' ',
                       as.character(ntree), ' ', as.character(p))
    old.stop <- nchar(old.text)

    post$treedraws$trees <- sub(old.text,
                                paste0(as.character(C*mc.cores*mc.ndpost),
                                       ' ', as.character(ntree), ' ',
                                       as.character(p)),
                                post$treedraws$trees)

    for(h in 1:C)
        for(i in 1:mc.cores)
            if(!(h==1 && i==1)) {

                post$yhat.shard[(i-1)*mc.ndpost+1:mc.ndpost, shard[[h]] ] <-
                    post.list[[h]][[i]]$yhat.train

                post$varcount <- rbind(post$varcount,
                                       post.list[[h]][[i]]$varcount)

                post$treedraws$trees <-
                    paste0(post$treedraws$trees,
                           substr(post.list[[h]][[i]]$treedraws$trees,
                                  old.stop+2,
                                  nchar(post.list[[h]][[i]]$treedraws$trees)))
            }

    attr(post, 'class') <- 'wbart'

    if(keeptrainfits) {
        post$x.train <- t(x.train)
        post$yhat.train <- predict(post,
                                   newdata=post$x.train, mc.cores=mc.cores)
    }

    if(length(x.test)>0) {
        post$x.test <- bartModelMatrix(x.test)
        post$yhat.test <-
            predict(post, newdata=post$x.test, mc.cores=mc.cores)
    }

    return(post)
}
