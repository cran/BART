
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

pbart=function(
x.train, y.train, x.test=matrix(0.0,0,0),
sparse=FALSE, theta=0, omega=1,
a=0.5, b=1, augment=FALSE, rho=NULL,
xinfo=matrix(0.0,0,0), usequants=FALSE,
cont=FALSE, rm.const=TRUE,
k=2.0, power=2.0, base=.95,
binaryOffset=NULL,
ntree=50L, numcut=100L,
ndpost=1000L, nskip=100L, keepevery=1L,
nkeeptrain=ndpost, nkeeptest=ndpost,
##nkeeptestmean=ndpost,
nkeeptreedraws=ndpost,
printevery=100L, transposed=FALSE
##treesaslists=FALSE
)
{
#--------------------------------------------------
#data
n = length(y.train)

if(length(binaryOffset)==0) binaryOffset=qnorm(mean(y.train))

if(!transposed) {
    temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                           cont=cont, xinfo=xinfo, rm.const=rm.const)
    x.train = t(temp$X)
    numcut = temp$numcut
    xinfo = temp$xinfo
    ## if(length(x.test)>0)
    ##         x.test = t(bartModelMatrix(x.test[ , temp$rm.const]))
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
if(length(rho)==0) rho <- p
if(length(rm.const)==0) rm.const <- 1:p
if(length(grp)==0) grp <- 1:p

#--------------------------------------------------
#set  nkeeps for thinning
if((nkeeptrain!=0) & ((ndpost %% nkeeptrain) != 0)) {
   nkeeptrain=ndpost
   cat('*****nkeeptrain set to ndpost\n')
}
if((nkeeptest!=0) & ((ndpost %% nkeeptest) != 0)) {
   nkeeptest=ndpost
   cat('*****nkeeptest set to ndpost\n')
}
## if((nkeeptestmean!=0) & ((ndpost %% nkeeptestmean) != 0)) {
##    nkeeptestmean=ndpost
##    cat('*****nkeeptestmean set to ndpost\n')
## }
if((nkeeptreedraws!=0) & ((ndpost %% nkeeptreedraws) != 0)) {
   nkeeptreedraws=ndpost
   cat('*****nkeeptreedraws set to ndpost\n')
}
#--------------------------------------------------
#prior
## nu=sigdf
## if(is.na(lambda)) {
##    if(is.na(sigest)) {
##       if(p < n) {
##          df = data.frame(t(x.train),y.train)
##          lmf = lm(y.train~.,df)
##          rm(df)
##          sigest = summary(lmf)$sigma
##       } else {
##          sigest = sd(y.train)
##       }
##    }
##    qchi = qchisq(1.0-sigquant,nu)
##    lambda = (sigest*sigest*qchi)/nu #lambda parameter for sigma prior
## }

## if(is.na(sigmaf)) {
##    tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree));
## } else {
##    tau = sigmaf/sqrt(ntree)
## }
#--------------------------------------------------
ptm <- proc.time()
#call
res = .Call("cpbart",
            n,  #number of observations in training data
            p,  #dimension of x
            np, #number of observations in test data
            x.train,   #p*n training data x
            y.train,   #n*1 training data y
            x.test,    #p*np test data x
            ntree,
            numcut,
            ndpost*keepevery,
            nskip,
            power,
            base,
            binaryOffset,
            3/(k*sqrt(ntree)),
            sparse,
            theta,
            omega,
            grp,
            a,
            b,
            rho,
            augment,
            nkeeptrain,
            nkeeptest,
            ##nkeeptestmean,
            nkeeptreedraws,
            printevery,
            ##treesaslists,
            xinfo
)
    
res$proc.time <- proc.time()-ptm
    
if(nkeeptrain>0) {
    ##res$yhat.train.mean <- NULL
    ##res$yhat.train.mean = res$yhat.train.mean+binaryOffset
    res$yhat.train = res$yhat.train+binaryOffset
    res$prob.train = pnorm(res$yhat.train)
    res$prob.train.mean <- apply(res$prob.train, 2, mean)
} else {
    res$yhat.train <- NULL
    ##res$yhat.train.mean <- NULL
}

if(np>0) {
    ##res$yhat.test.mean <- NULL
    ##res$yhat.test.mean = res$yhat.test.mean+binaryOffset
    res$yhat.test = res$yhat.test+binaryOffset
    res$prob.test = pnorm(res$yhat.test)
    res$prob.test.mean <- apply(res$prob.test, 2, mean)
} else {
    res$yhat.test <- NULL
    ##res$yhat.test.mean <- NULL
}

if(nkeeptreedraws>0) ## & !treesaslists)
    names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]

dimnames(res$varcount)[[2]] = as.list(dimnames(x.train)[[1]])
dimnames(res$varprob)[[2]] = as.list(dimnames(x.train)[[1]])
res$varcount.mean <- apply(res$varcount, 2, mean)
res$varprob.mean <- apply(res$varprob, 2, mean)
res$rm.const <- rm.const
res$binaryOffset=binaryOffset
attr(res, 'class') <- 'pbart'
return(res)
}
