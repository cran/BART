
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

gewekediag <-  function (x, frac1 = 0.1, frac2 = 0.5)
{
    if (frac1 < 0 || frac1 > 1) stop("frac1 invalid")

    if (frac2 < 0 || frac2 > 1) stop("frac2 invalid")

    if (frac1 + frac2 > 1) stop("start and end sequences are overlapping")

    end. <- nrow(x)
    xstart <- c(1, floor(end. - frac2 * (end. - 1)))
    xend <- c(ceiling(1 + frac1 * (end. - 1)), end.)
    y.variance <- y.mean <- vector("list", 2)

    for (i in 1:2) {
        y <- x[xstart[i]:xend[i], ]
        y.mean[[i]] <- apply(y, 2, mean)
        y.variance[[i]] <- spectrum0ar(y)$spec/(xend[i]-xstart[i]+1)
    }

    z <- (y.mean[[1]] - y.mean[[2]])/sqrt(y.variance[[1]] + y.variance[[2]])
    out <- list(z = z, frac = c(frac1, frac2))

    return(out)
}
