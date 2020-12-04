
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

gbart=function(
               x.train, y.train,
               x.test=matrix(0,0,0), type='wbart',
               ntype=as.integer(
                   factor(type, levels=c('wbart', 'pbart', 'lbart'))),
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
               mc.cores = 1L, nice = 19L, seed = 99L
               )
{
    if(is.na(ntype))
        stop("type argument must be set to either 'wbart', 'pbart' or 'lbart'")

    n = length(y.train)

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
        grp <- temp$grp
        rm(temp)
    }
    else {
        rm.const <- NULL
        grp <- NULL
    }

    if(n!=ncol(x.train))
        stop('The length of y.train and the number of rows in x.train must be identical')

    p = nrow(x.train)
    np = ncol(x.test)
    if(length(rho)==0) rho=p
    if(length(rm.const)==0) rm.const <- 1:p
    if(length(grp)==0) grp <- 1:p

    check <- unique(sort(y.train))

    if(length(check)==2) {
        if(!all(check==0:1))
            stop('Binary y.train must be coded as 0 and 1')
        if(type=='wbart')
            stop("The outcome is binary so set type to 'pbart' or 'lbart'")
    }

    ## check <- c('wbart', 'pbart', 'lbart')

    ## if(!(type %in% check))
    ##     stop("type argument must be set to either 'wbart', 'pbart' or 'lbart'")

    if(length(offset)==0) {
        offset=mean(y.train)
        if(type=='pbart') offset=qnorm(offset)
        else if(type=='lbart') offset=qlogis(offset)
    }

    if(type=='wbart') {
        y.train = y.train-offset

        if(!is.na(sigest) && !is.na(lambda) && lambda==0) {
            ##no op: sigma is fixed and known at given sigest value
        }
        else if(is.na(lambda)) {
            if(is.na(sigest)) {
                if(p < n)
                    sigest = summary(lm(y.train~.,
                                        data.frame(t(x.train),y.train)))$sigma
                else sigest = sd(y.train)
            }
            qchi = qchisq(1-sigquant, sigdf)
            lambda = (sigest^2)*qchi/sigdf #lambda parameter for sigma prior
        } else {
            sigest=sqrt(lambda)
        }

        if(is.na(tau.num)) {
            tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
        } else {
            tau=tau.num/(k*sqrt(ntree))
        }
    } else {
        lambda=1
        sigest=1
        tau=tau.num/(k*sqrt(ntree))
        ## tau=1-tau.interval

        ## if(type=='pbart')
        ##     tau=qnorm(1-0.5*tau)/(k*sqrt(ntree))
        ## else if(type=='lbart')
        ##     tau=qlogis(1-0.5*tau)/(k*sqrt(ntree))
    }

    ## hot deck missing imputation
    ## must be conducted here since it would
    ## cause trouble with multi-threading on the C++ side

    check=(np>0 && np==n)

    for(i in 1:n)
        for(j in 1:p) {
            if(check) check=((is.na(x.train[j, i]) && is.na(x.test[j, i])) ||
                             (!is.na(x.train[j, i]) && !is.na(x.test[j, i]) &&
                              x.train[j, i]==x.test[j, i]))

            while(is.na(x.train[j, i])) {
                h=sample.int(n, 1)
                x.train[j, i]=x.train[j, h]
            }
        }

    if(check) x.test=x.train
    else if(np>0) {
        for(i in 1:np)
            for(j in 1:p)
                while(is.na(x.test[j, i])) {
                    h=sample.int(np, 1)
                    x.test[j, i]=x.test[j, h]
                }
    }

    ## if(hotdeck) ## warnings are suppressed with mc.gbart anyways
    ##     warning('missing elements of x imputed with hot decking')

    if(.Platform$OS.type!='unix') hostname <- FALSE
    else if(hostname)
        hostname <- system('hostname', intern=TRUE)

    ptm <- proc.time()

    res = .Call("cgbart",
                ntype, ##as.integer(factor(type, levels=check))-1,
                n,  #number of observations in training data
                p,  #dimension of x
                np, #number of observations in test data
                x.train,   #pxn training data x
                y.train,   #pxn training data x
                x.test,    #p*np test data x
                ntree,
                numcut,
                ndpost*keepevery,
                nskip,
                keepevery,
                power,
                base,
                offset,
                tau,
                sigdf,
                lambda,
                sigest,
                w,
                sparse,
                theta,
                omega,
                grp,
                a,
                b,
                rho,
                augment,
                printevery,
                xinfo
                )

    res$proc.time <- proc.time()-ptm
    res$hostname <- hostname

    Y=t(matrix(y.train, nrow=n, ncol=ndpost))

    if(type=='wbart') {
        res$yhat.train.mean <- apply(res$yhat.train, 2, mean)
        SD=matrix(res$sigma[-(1:nskip)], nrow=ndpost, ncol=n)
        ##CPO=1/apply(1/dnorm(Y, res$yhat.train, SD), 2, mean)
        log.pdf=dnorm(Y, res$yhat.train, SD, TRUE)
        res$sigma.mean=mean(SD[ , 1])
    }
    else {
        if(type=='pbart') res$prob.train = pnorm(res$yhat.train)
        else if(type=='lbart') res$prob.train = plogis(res$yhat.train)

        ##CPO=1/apply(1/dbinom(Y, 1, res$prob.train), 2, mean)
        log.pdf=dbinom(Y, 1, res$prob.train, TRUE)

        res$prob.train.mean <- apply(res$prob.train, 2, mean)
    }

    min.log.pdf=t(matrix(apply(log.pdf, 2, min), nrow=n, ncol=ndpost))
    log.CPO=log(ndpost)+min.log.pdf[1, ]-
        log(apply(exp(min.log.pdf-log.pdf), 2, sum))
    res$LPML=sum(log.CPO)
    ##res$CPO=exp(log.CPO)
    ##res$LPML=sum(log(CPO))

    if(np>0) {
        if(type=='wbart')
            res$yhat.test.mean <- apply(res$yhat.test, 2, mean)
        else {
            if(type=='pbart') res$prob.test = pnorm(res$yhat.test)
            else if(type=='lbart') res$prob.test = plogis(res$yhat.test)

            res$prob.test.mean <- apply(res$prob.test, 2, mean)
        }
    }

    res$ndpost = ndpost
    res$offset = offset
    names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]
    dimnames(res$varcount)[[2]] = as.list(dimnames(x.train)[[1]])
    dimnames(res$varprob)[[2]] = as.list(dimnames(x.train)[[1]])
    res$varcount.mean <- apply(res$varcount, 2, mean)
    res$varprob.mean <- apply(res$varprob, 2, mean)
    res$rm.const <- rm.const
    attr(res, 'class') <- type
    return(res)
}
