
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

spectrum0ar <- function(x)
{
  v0 <- order <- numeric(ncol(x))
  names(v0) <- names(order) <- colnames(x)
  z <- 1:nrow(x)

  for (i in 1:ncol(x)) {
      lm.out <- lm(x[,i] ~ z)

      if (identical(all.equal(sd(residuals(lm.out)), 0), TRUE)) {
          v0[i] <- 0
          order[i] <- 0
      }
      else {
          ar.out <- ar(x[,i], aic=TRUE)
          v0[i] <- ar.out$var.pred/(1 - sum(ar.out$ar))^2
          order[i] <- ar.out$order
      }
  }

  return(list(spec=v0, order=order))
}
