
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

mc.wbart.gse <- function(
                    x.train, y.train, 
                    P=50L, ## number of permutations 
                    R=5L,  ## number of replicates 
                    ntree=20L, numcut=100L,  
                    C=1, alpha=0.05,
                    k=2.0, power=2.0, base=0.95,
                    ndpost=2000L, nskip=100L,
                    printevery=100L, keepevery=1L, keeptrainfits=FALSE,
                    seed = 99L,    
                    mc.cores = 2L, 
                    nice=19L       
                    )
{
    if(.Platform$OS.type!='unix')
        stop('parallel::mcparallel/mccollect do not exist on windows')

    set.seed(seed)

    N <- length(y.train)
    K <- ncol(x.train)
    
    perm <- matrix(runif(N*P), nrow=N, ncol=P)
    
    prob <- matrix(nrow=P, ncol=K)
    
    h <- 1
    
    for(i in 0:P) {
        if(i==0) y. <- y.train
        else y. <- y.train
        
        tmp2 <- matrix(nrow=R, ncol=K)
        
        for(j in 1:R) {
            tmp1 <- mc.wbart(x.train=x.train, y.train=y.,
                                 k=k, power=power, base=base,
                                 ntree=ntree, numcut=numcut,
                                 ndpost=ndpost, nskip=nskip,
                                 printevery=printevery,
                                 keepevery=keepevery,
                                 keeptrainfits=keeptrainfits,
                                 seed=h, mc.cores=mc.cores, nice=nice)$varcount

            tmp2[j, ] <- apply(tmp1, 2, mean)

            h <- h+1
        }

        tmp1 <- apply(tmp2, 2, mean)
        tmp1 <- tmp1/sum(tmp1)
        
        if(i==0) varcount <- tmp1
        else prob[i, ] <- tmp1
    }

    mu.k <- apply(prob, 2, mean)
    sd.k <- apply(prob, 2, sd)

    cov.prob <- double(K)

    iter <- 0

    while(min(cov.prob)<(1-alpha)) {
        if(iter>0) {
            C <- C*1.01
            cov.prob <- cov.prob*0
        }
        
        for(i in 1:P) for(j in 1:K) {
            cov.prob[j] <- cov.prob[j]+(prob[i, j]<=(mu.k[j]+C*sd.k[j]))/P
        }

        iter <- iter+1
    }

    if(iter==1) 
        warning('Algorithm stopped at iteration 1.  Try again with a smaller C.')
    
    return(list(which=which(varcount>(mu.k+C*sd.k)), prob=varcount,
                C=C, mu.k=mu.k, sd.k=sd.k, iter=iter, perm.prob=prob))
}
               
