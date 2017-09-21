\name{recur.pre.bart}
\title{Data construction for recurrent events with BART}
\alias{recur.pre.bart}
\description{
Recurrent event data contained in \eqn{(t_1,\delta_1, ..., t_k,\delta_k, x)} must be translated to data
suitable for the BART model; see \code{recur.bart} for more details.
}
\usage{
recur.pre.bart( times, delta, x.train=NULL, probs=c(0.15, 0.15), baseline=FALSE )
}
\arguments{
   \item{times}{
     Matrix of time to event or right-censoring.\cr
   }

   \item{delta}{
     Matrix of event indicators: 1 is an event while 0 is censored.\cr
   }

   \item{x.train}{
   Explanatory variables for training (in sample) data.\cr
   If provided, must be a matrix
   with (as usual) rows corresponding to observations and columns to variables.\cr
   }

   \item{probs}{ \code{tx.test} is generated from \code{x.train} in a
     similar fashion to \code{tx.train}.  However, observations beyond
     censoring contain \code{NA} for columns two and three containing
     \eqn{v(t)} and \eqn{N(t-)} respectively.  In order to preform
     prediction with \code{tx.test}, then these missing values need to be
     handled accordingly.  We provide an automated method for this which
     is returned as \code{X.base}.  It would be difficult to provide an
     efficient method that randomly generates \eqn{v(t)} and \eqn{N(t-)}
     since they are structurally dependent.  Rather, we
     deterministically generate them in such a fashion that they define
     a "middle" baseline, i.e. the median \eqn{N(t-)} as a function of
     \eqn{t} for non-censored subjects (calculated from \code{tx.train});
     the median is adjusted upward as the censoring increases with time
     (if it does at all) to a maximum value controlled by the second
     argument contained in \code{probs} since otherwise the "low" risk
     subjects would be over-represented.  Each subject's previous value
     is carried forward unless their computed sojourn time (which is a
     function of the generated \eqn{N(t-)}) is longer than the the
     median \eqn{v(t)} as a function of \eqn{t} (which is also
     calculated from \code{tx.train}); the median is adjusted downward as
     the censoring increases with time (if it does at all) to a maximum
     value controlled by the first argument contained in \code{probs}. }

   \item{baseline}{experimental feature: still needs work}
}

\value{
   \code{recur.pre.bart} returns a list.
Besides the items listed below, the list has
 a \code{times} component giving the unique times and \code{K} which is the number of
 unique times.
   \item{y.train}{A vector of binary responses.}
   \item{tx.train}{A matrix with the rows of the training data.}
   \item{tx.test}{Generated from \code{x.train} (see discussion above included in
     the argument \code{probs}).}
   \item{X.base}{Generated from \code{x.train} (see discussion above included in
     the argument \code{probs}).}
   \item{sojourn}{Generated from \code{x.train} (see discussion above included in
     the argument \code{probs}).}
   \item{pattern}{Generated from \code{x.train} (see discussion above included in
     the argument \code{probs}).}
}
\references{
Sparapani, Rein, Tarima, Jackson, Meurer (2017).  Nonparametric recurrent
events analysis with BART and an application to the hospital admissions
of patients with diabetes.  MCW Biostatistics Technical Report 64
\url{http://www.mcw.edu/FileLibrary/Groups/Biostatistics/TechReports/TechReports5175/tr064.pdf}
}
\author{
Rodney Sparapani: \email{rsparapa@mcw.edu}
}
\seealso{
\code{\link{recur.bart}}
}
\examples{

data(bladder)
subset <- -which(bladder1$stop==0)
bladder0 <- bladder1[subset, ]
id <- unique(sort(bladder0$id))
N <- length(id)
K <- max(bladder0$enum)

times <- matrix(0, nrow=N, ncol=K)
dimnames(times)[[1]] <- paste0(id)

delta <- matrix(0, nrow=N, ncol=K)
dimnames(delta)[[1]] <- paste0(id)

x.train <- matrix(NA, nrow=N, ncol=2*K+1) ## add time-dependent cols too
dimnames(x.train)[[1]] <- paste0(id)
dimnames(x.train)[[2]] <- c('rx', rep(c('number', 'size'), K))

for(i in 1:N) {
    h <- id[i]

    for(j in 1:K) {
        k <- which(bladder0$id==h & bladder0$enum==j)

        if(length(k)==1) {
            times[i, j] <- bladder0$stop[k]
            delta[i, j] <- (bladder0$status[k]==1)*1

            if(j==1) {
                x.train[i, 1] <- as.numeric(bladder0$treatment[k])
                x.train[i, 2] <- bladder0$number[k]
                x.train[i, 3] <- bladder0$size[k]
            }
            else if(delta[i, j]==1) {
                if(bladder0$rtumor[k]!='.')
                    x.train[i, 2*j] <- as.numeric(bladder0$rtumor[k])
                if(bladder0$rsize[k]!='.')
                    x.train[i, 2*j+1] <- as.numeric(bladder0$rsize[k])
            }
        }
    }
}

\dontrun{
pre <- recur.pre.bart(times=times, delta=delta, x.train=x.train)

J <- nrow(pre$tx.train)

for(j in 1:J) {
    if(pre$tx.train[j, 3]>0) {
        pre$tx.train[j, 5] <- pre$tx.train[j, 5+pre$tx.train[j, 3]*2]
        pre$tx.train[j, 6] <- pre$tx.train[j, 6+pre$tx.train[j, 3]*2]
    }
}

pre$tx.train <- pre$tx.train[ , 1:6]

## tx.test is a little more work since you also need to set v and N
J <- nrow(pre$tx.test)

for(j in 1:J) {
    if(is.na(pre$tx.test[j, 3])) {
        pre$tx.test[j, 2] <- pre$tx.test[j, 1]-t.
        pre$tx.test[j, 3] <- N.
    }
    else if(pre$tx.test[j, 1]==1) {
        t. <- 0
        N. <- 0
    }
    else if(pre$tx.test[j, 3]>N.) {
        t. <- pre$tx.test[j, 1]
        N. <- pre$tx.test[j, 3]
        pre$tx.test[j, 5] <- pre$tx.test[j, 5+pre$tx.test[j, 3]*2]
        pre$tx.test[j, 6] <- pre$tx.test[j, 6+pre$tx.test[j, 3]*2]
    }
}

pre$tx.test <- pre$tx.test[ , 1:6]

## since BART is essentially a nonparametric function, handling
## missing data can, in certain circumstances, be as easy as creating
## a missing "category" for continuous and ordinal variables
## in the bladder1 ex, both number and size are recorded as integers
## from 1 to 8 so a reasonable missing category might be 0 or 9
## however, we really want to create a generous amount of
## separation between observed and missing data to give BART a
## reasonable chance of splitting between them; but, perhaps, not too
## generous otherwise you may sacrfice performance within the
## observed data range; the usual caveats apply and it is always
## a good idea to explore more sophisticated methods (see below)    
pre$tx.train[which(is.na(pre$tx.train[ , 5])), 5] <- -2
pre$tx.train[which(is.na(pre$tx.train[ , 6])), 6] <- -2
pre$tx.test[which(is.na(pre$tx.test[ , 5])), 5] <- -2
pre$tx.test[which(is.na(pre$tx.test[ , 6])), 6] <- -2
    
## or impute missing data via your favorite method
## such as Sequential BART (see ref above)
## library(MASS)
## library(Rlab)
## library(msm)
## source('~/sequential_bart/load-pbart.R')
## set.seed(21)
## train <- BartMI(bartmpicstem, xx=pre$tx.train, yy=NULL, datatype=rep(0, 6),
##                type=0, numskip=20, burn=1000)
## coarsen the imputed data same way as observed example data
## train$imputed5[which(train$imputed5[ , 5]<1), 5] <- 1
## train$imputed5[which(train$imputed5[ , 5]>8), 5] <- 8
## train$imputed5[ , 5] <- round(train$imputed5[ , 5])
## train$imputed5[which(train$imputed5[ , 6]<1), 6] <- 1
## train$imputed5[which(train$imputed5[ , 6]>8), 6] <- 8
## train$imputed5[ , 6] <- round(train$imputed5[ , 6])
## pre$tx.train <- train$imputed5
## dimnames(pre$tx.train)[[2]] <- c('t', 'v', 'N', 'rx', 'number', 'size')

## for Friedman's partial dependence, we need to estimate the whole cohort
## at each treatment assignment (and, average over those)
## Pl treatment=1, B6 treatment=2, Th treatment=3
pre$tx.test <- rbind(pre$tx.test, pre$tx.test, pre$tx.test)
pre$tx.test[ , 4] <- rep(1:3, each=J)
    
set.seed(99)
post <- recur.bart(y.train=pre$y.train, x.train=pre$tx.train, x.test=pre$tx.test)
## depending on your performance, you may want to run in parallel if available
## post <- mc.recur.bart(y.train=pre$y.train, x.train=pre$tx.train, 
##                       x.test=pre$tx.test, mc.cores=5, seed=99)

K <- post$K
M <- 1000
haz.Pl <- matrix(0, nrow=M, ncol=K)
haz.B6 <- matrix(0, nrow=M, ncol=K)
haz.Th <- matrix(0, nrow=M, ncol=K)

for(i in 1:N) {
    j <- (i-1)*K+1:K
    haz.Pl <- haz.Pl+post$haz.test[1:M, j]/N
    haz.B6 <- haz.B6+post$haz.test[1:M, j+N*K]/N
    haz.Th <- haz.Th+post$haz.test[1:M, j+N*K*2]/N
}

haz.Pl.mu <- apply(haz.Pl, 2, mean)
haz.B6.mu <- apply(haz.B6, 2, mean)
haz.Th.mu <- apply(haz.Th, 2, mean)

plot(post$times, haz.Pl.mu, type='l',
     main='bladder1 ex: Placebo (black), B6 (red), Thiotepa (blue)',
     ylab=expression(lambda(t)), xlab='t (months)',
     ylim=c(0, max(c(haz.Pl.mu, haz.B6.mu, haz.Th.mu))))
points(post$times, haz.B6.mu, type='l', col='red')
points(post$times, haz.Th.mu, type='l', col='blue')

plot(post$times, apply(1*((haz.Th/haz.Pl) < 1), 2, mean),
     type='l', ylim=c(0, 1),
     main='bladder1 ex: Thiotepa vs. Placebo',
     ylab=expression(P(RR(t)<1)), xlab='t (months)')
     points(c(0, max(post$times)), c(0.5, 0.5), type='l', col='gray')

rr.Th.Pl.mu <- apply(haz.Th/haz.Pl, 2, mean)
rr.Th.Pl.lo <- apply(haz.Th/haz.Pl, 2, quantile, probs=0.25)
rr.Th.Pl.up <- apply(haz.Th/haz.Pl, 2, quantile, probs=0.75)
plot(post$times, rr.Th.Pl.mu, type='l', log='y',
     main='bladder1 ex: Thiotepa vs. Placebo with 0.5 CI',
     ylab=expression(RR(t)), xlab='t (months)', ylim=c(0.625, 1.6))
points(post$times, rr.Th.Pl.lo, type='l', col='gray')
points(post$times, rr.Th.Pl.up, type='l', col='gray')
points(c(0, max(post$times)), c(1, 1), type='l')

rr.Th.B6.mu <- apply(haz.Th/haz.B6, 2, mean)
rr.Th.B6.lo <- apply(haz.Th/haz.B6, 2, quantile, probs=0.25)
rr.Th.B6.up <- apply(haz.Th/haz.B6, 2, quantile, probs=0.75)
plot(post$times, rr.Th.B6.mu, type='l', log='y',
     main='bladder1 ex: Thiotepa vs. Vitamin B6 with 0.5 CI',
     ylab=expression(RR(t)), xlab='t (months)', ylim=c(0.625, 1.6))
points(post$times, rr.Th.B6.lo, type='l', col='gray')
points(post$times, rr.Th.B6.up, type='l', col='gray')
points(c(0, max(post$times)), c(1, 1), type='l')

plot(post$times, apply(1*((haz.Th/haz.B6) < 1), 2, mean),
     type='l', ylim=c(0, 1),
     main='bladder1 ex: Thiotepa vs. Vitamin B6',
     ylab=expression(P(RR(t)<1)), xlab='t (months)')
     points(c(0, max(post$times)), c(0.5, 0.5), type='l', col='gray')


}
}
\keyword{data construction}
