
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


recur.pre.bart <- function(
                           times,
                           ## matrix of recur times (start times)

                           delta,
                           ## matrix of indicators: 0=no event, 1=event

                           x.train=NULL,
                           ## matrix of covariate regressors
                           ## can be NULL, i.e. KM analog

                           tstop=NULL,
                           ## for non-instantaneous events, this the
                           ## matrix of event stop times, i.e., between
                           ## times[i, j] and tstop[i, j] subject i is
                           ## not in the risk set for a recurrent event
                           ## N.B. NOT for counting process notation

                           last.value=TRUE
                           ) {
    ## currently does not handle time dependent Xs
    ## can be extended later
    ## most likely via the alternative counting process notation

    N <- nrow(times)

    if(N!=nrow(delta))
        stop('The number of rows in times and delta must be identical')

    if(length(x.train)>0 && N!=nrow(x.train))
        stop('The number of rows in times and x.train, if any, must be identical')

    J <- ncol(times)

    times <- cbind(times, apply(times, 1, max, na.rm=TRUE))

    dimnames(times)[[2]][J+1] <- 'stop'

    events <- unique(sort(c(times, tstop)))
    ## time grid of events including stop times

    if(events[1]==0) events <- events[-1]

    K <- length(events)

    y.train <- integer(N) ## y.train is at least N long

    k <- 1
    id <- 0
    for(i in 1:N)
        for(j in 1:K)
            if(events[j] <= times[i, J+1]) {
                y.train[k] <- 0
                id[k] <- i

                for(h in 1:J) if(!is.na(times[i, h]) && !is.na(delta[i, h]) &&
                                 y.train[k]==0 && events[j]==times[i, h] && delta[i, h]==1)
                                  y.train[k] <- 1

                k <- k+1
            }

    m <- length(y.train)

    ##binaryOffset <- qnorm(mean(y.train))
    
    if(length(x.train)==0) {
        p <- 0
        n <- 1

        X.train <- matrix(nrow=m, ncol=3)

        dimnames(X.train)[[2]] <- c('t', 'v', 'N')
    } else {
        p <- ncol(x.train)

        ##if(length(x.test)>0) n <- nrow(x.test)

        X.train <- matrix(nrow=m, ncol=p+3)

        if(length(dimnames(x.train)[[2]])>0)
            dimnames(X.train)[[2]] <- c('t', 'v', 'N', dimnames(x.train)[[2]])
        else dimnames(X.train)[[2]] <- c('t', 'v', 'N', paste0('x', 1:p))
    }

    X.test <- matrix(nrow=N*K, ncol=ncol(X.train))

    dimnames(X.test)[[2]] <- dimnames(X.train)[[2]]

    k <- 1

    for(i in 1:N) {
        n.t <- 0
        t.0 <- 0

        for(j in 1:K) if(events[j] <= times[i, J+1]) {
                          X.train[k, 1:3] <- c(events[j], events[j]-t.0, n.t)

                          for(h in 1:J) if(!is.na(times[i, h]) && !is.na(delta[i, h]) &&
                                           events[j]==times[i, h] && delta[i, h]==1) {
                                            n.t <- n.t+1
                                            t.0 <- events[j]
                                        }

                          if(p>0) X.train[k, 4:(3+p)] <- x.train[i, ]

                          k <- k+1
                      }
    }

    ## generate X.test from X.train with v(t) & N(t-) as NA beyond follow-up
    k <- 1

    for(i in 1:N) {
        n.t <- 0
        t.0 <- 0

        for(j in 1:K)  {
            X.test[k, 1:3] <- c(events[j], events[j]-t.0, n.t)

            for(h in 1:J) {
                if(!is.na(times[i, h]) && !is.na(delta[i, h]) &&
                   events[j]==times[i, h] && delta[i, h]==1) {
                    n.t <- n.t+1
                    t.0 <- events[j]
                }
                else if(!last.value & events[j] == times[i, J+1]) {
                    n.t <- NA
                    t.0 <- NA
                }
            }

            if(p>0) X.test[k, 4:(3+p)] <- x.train[i, ]

            k <- k+1
        }
    }

    ## remove training entries not in the risk-set
    h <- TRUE
    if(length(tstop)>0) {
        for(k in 1:m) {
            i <- id[k]
            t <- X.train[k, 1]
            j <- 1
            h[k] <- TRUE
            while(h[k] & j<=J) {
                h[k] <- (tstop[i, j]==0 | !(times[i, j]<t & t<=tstop[i, j]))
                j <- j+1
            }
        }

        y.train <- y.train[h]
        X.train <- X.train[h, ]
    }

    return(list(y.train=y.train, tx.train=X.train, tx.test=X.test,
                times=events, K=K)) ##, binaryOffset=binaryOffset))
}
