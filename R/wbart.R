
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

wbart=function( 
x.train, y.train, x.test=matrix(0.0,0,0),
sigest=NA, sigdf=3, sigquant=.90,
k=2.0,
power=2.0, base=.95,
sigmaf=NA,
lambda=NA,
fmean=mean(y.train),
w=rep(1,length(y.train)),
ntree=200, numcut=100,
ndpost=1000, nskip=100,
nkeeptrain=ndpost,nkeeptest=ndpost,
nkeeptestmean=ndpost,
nkeeptreedraws=ndpost,
printevery=100, transposed=FALSE,
treesaslists=FALSE
)
{
#--------------------------------------------------
nd = ndpost
burn = nskip
#--------------------------------------------------
#data 
n = length(y.train)

if(!transposed) {
    x.train = t(x.train)
    x.test = t(x.test)
}

if(n!=ncol(x.train))
    stop('The length of y.train and the number of rows in x.train must be identical')
    
p = nrow(x.train)
np = ncol(x.test)

y.train = y.train-fmean
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
if((nkeeptestmean!=0) & ((ndpost %% nkeeptestmean) != 0)) {
   nkeeptestmean=ndpost
   cat('*****nkeeptestmean set to ndpost\n')
}
if((nkeeptreedraws!=0) & ((ndpost %% nkeeptreedraws) != 0)) {
   nkeeptreedraws=ndpost
   cat('*****nkeeptreedraws set to ndpost\n')
}
#--------------------------------------------------
#prior
nu=sigdf
if(is.na(lambda)) {
   if(is.na(sigest)) {
      if(p < n) {
         df = data.frame(t(x.train),y.train)
         lmf = lm(y.train~.,df)
         sigest = summary(lmf)$sigma
      } else {
         sigest = sd(y.train)
      }
   }
   qchi = qchisq(1.0-sigquant,nu)
   lambda = (sigest*sigest*qchi)/nu #lambda parameter for sigma prior
} else {
   sigest=sqrt(lambda)
}

if(is.na(sigmaf)) {
   tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree));
} else {
   tau = sigmaf/sqrt(ntree)
}
#--------------------------------------------------
#call
res = .Call("cwbart",
            n,  #number of observations in training data
            p,  #dimension of x
            np, #number of observations in test data
            x.train,   #pxn training data x
            y.train,   #pxn training data x
            x.test,   #p*np test data x
            ntree,
            numcut,
            nd,
            burn,
            power,
            base,
            tau,
            nu,
            lambda,
            sigest,
            w,
            nkeeptrain,
            nkeeptest,
            nkeeptestmean,
            nkeeptreedraws,
            printevery,
            treesaslists
)
res$yhat.train.mean = res$yhat.train.mean+fmean
res$yhat.train = res$yhat.train+fmean
res$yhat.test.mean = res$yhat.test.mean+fmean
res$yhat.test = res$yhat.test+fmean
##res$nkeeptreedraws=nkeeptreedraws
return(res)
}
