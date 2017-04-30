
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


## you call this function before bart()
## this function takes traditional time/delta
## recurrent event variables and regressors (if any)
## and it constructs the corresponding
## tx.train, y.train and tx.test appropriate
## for use with pbart()

recur.pre.bart <- function(
                      times,
                      ## matrix of recur times

                      delta,
                      ## matrix of indicators: 0=no event, 1=event

                      x.train=NULL,
                      ## matrix of covariate regressors
                      ## can be NULL, i.e. KM analog

                      probs=c(0.15, 0.15)
                      ## the maximum amount to adjust the quantile
                      ## for the middle pattern due to censoring
                      ## first decreases v(t) and the second increases N(t-)
, baseline=FALSE
                      ##x.test=NULL
                      ## matrix of covariate regressors at tx.test settings
                      ## does nothing for now since there is no obvious basis for v(t) and N(t-)
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

    times <- cbind(times, apply(times, 1, max))

    dimnames(times)[[2]][J+1] <- 'stop'

    events <- unique(sort(times))
    ## time grid of events including censoring times

    if(events[1]==0) events <- events[-1]

    K <- length(events)

    y.train <- integer(N) ## y.train is at least N long

    k <- 1

    for(i in 1:N) for(j in 1:K) if(events[j] <= times[i, J+1]) {
        y.train[k] <- 0

        for(h in 1:J) if(y.train[k]==0 & events[j]==times[i, h] & delta[i, h]==1)
            y.train[k] <- 1

        k <- k+1
    }

    m <- length(y.train)

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

            for(h in 1:J) if(events[j]==times[i, h] & delta[i, h]==1) {
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
                if(events[j]==times[i, h] & delta[i, h]==1) {
                    n.t <- n.t+1
                    t.0 <- events[j]
                }
                else if(events[j] == times[i, J+1]) {
                    n.t <- NA
                    t.0 <- NA
                }
            }

            if(p>0) X.test[k, 4:(3+p)] <- x.train[i, ]

            k <- k+1
        }
    }

    ## generate X.base from X.train with v(t) & N(t-): baseline of the "middle" pattern
if(baseline) {
    X.base <- cbind(X.test)

    dimnames(X.base)[[2]] <- dimnames(X.train)[[2]]

    sojourn <- double(K)
    pattern <- double(K)

    for(j in 1:K) {
        h <- seq(j, N*K, K)
        sojourn[j] <- quantile(X.base[h, 2], na.rm=TRUE,
                               probs=0.5-probs[1]*(1-mean(1*(!is.na(X.base[h, 2])))))
        pattern[j] <- round(quantile(X.base[h, 3], na.rm=TRUE,
                                     probs=0.5+probs[2]*(1-mean(1*(!is.na(X.base[h, 3]))))))

        if(j>1) {
            if(pattern[j-1]<pattern[j]) pattern[j] <- pattern[j-1]+1
            else if(pattern[j-1]>pattern[j]) pattern[j] <- pattern[j-1]
        }
    }

    for(i in 1:N) {
        t.0 <- 0

        for(j in 1:K) {
            h <- (i-1)*K+j
            if(is.na(X.base[h, 3])) {
                if(X.base[h-1, 3]>pattern[j]) X.base[h, 3] <- X.base[h-1, 3]
                else if(X.base[h-1, 3]<pattern[j] & (X.base[h-1, 1]-t.0)>=sojourn[j]) {
                ##else if(X.base[h-1, 3]<pattern[j]) {
                    t.0 <- X.base[h-1, 1]
                    X.base[h, 3] <- X.base[h-1, 3]+1
                }
                else X.base[h, 3] <- X.base[h-1, 3]

                X.base[h, 2] <- X.base[h, 1]-t.0
            }
            else if(j>1) if(X.base[h, 3]>X.base[h-1, 3]) t.0 <- X.base[h-1, 1]
        }
    }
    } else {
        X.base <- NULL
        pattern <- NULL
        sojourn <- NULL
        }

## automated X.test creation is not feasible since there is no obvious basis for v(t) and N(t-)
    ## if(p==0 | length(x.test)>0) {
    ##     X.test <- matrix(nrow=K*n, ncol=p+3, dimnames=dimnames(X.train))

    ##     ## we can summarize this, but this is only one of many possibilities
    ##     for(i in 1:n) for(j in 1:K) {
    ##         X.test[(i-1)*K+j, 1:(J+1)] <- c(rep(events[j], J), 0)
    ##         if(p>0) X.test[(i-1)*K+j, (J+2):(J+1+p)] <- x.test[i, ]
    ##     }
    ## }
    ## else X.test <- matrix(0.0, 0L, 0L)

    return(list(y.train=y.train, tx.train=X.train, tx.test=X.test, tx.base=X.base,
                times=events, K=K, pattern=pattern, sojourn=sojourn))

    ##return(list(y.train=y.train, X.train=X.train, X.test=matrix(0.0, 0L, 0L), times=events, K=K))
                ##X.test=data.matrix(X.test), times=events, K=K))
}
