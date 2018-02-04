
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017 Rodney Sparapani

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

stratrs <- function(y, C=floor(length(y)/2000))
{
    strat <- unique(sort(y))
    rs <- integer(length(y))
    
    for(i in strat) {
        N <- sum(y==i)
        rs[y==i] <- sample(1:N, N)%%C
    }
    
    return(rs+1)
}

## set.seed(12)
## x <- rbinom(25000, 1, 0.1)
## a <- stratrs(x)
## table(a, x)
## z <- pmin(rpois(25000, 0.8), 5)
## b <- stratrs(z)
## table(b, z)
