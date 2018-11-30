
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017-2018 Robert McCulloch and Rodney Sparapani

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

pwbart = function(
   x.test,		#x matrix to predict at
   treedraws,		#$treedraws from wbart
   mu=0,		#mean to add on
   mc.cores=1L,         #thread count
   transposed=FALSE,	
   dodraws=TRUE,
   nice=19L             #mc.pwbart only	
)
{
if(!transposed) x.test <- t(bartModelMatrix(x.test))

p <- length(treedraws$cutpoints)

if(p!=nrow(x.test))
    stop(paste0('The number of columns in x.test must be equal to ', p))

res = .Call("cpwbart",
   treedraws,	#trees list
   x.test,      #the test x
   mc.cores   	#thread count
)
if(dodraws) return(res$yhat.test+mu)
else return(apply(res$yhat.test, 2, mean)+mu)
}
