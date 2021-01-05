\name{recur.pre.bart}
\title{Data construction for recurrent events with BART}
\alias{recur.pre.bart}
\description{
Recurrent event data contained in \eqn{(t_1,\delta_1, ..., t_k,\delta_k, x)} must be translated to data
suitable for the BART model; see \code{recur.bart} for more details.
}
\usage{
recur.pre.bart( times, delta, x.train=NULL, tstop=NULL, last.value=TRUE )
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

   \item{tstop}{ For non-instantaneous events, this the matrix of event
     stop times, i.e., between \code{times[i, j]} and \code{tstop[i, j]}
     subject \code{i} is not in the risk set for a recurrent event.
     N.B. This is NOT for counting process notation. \cr }

   \item{last.value}{If \code{last.value=TRUE}, then the sojourn time,
     \code{v}, and the number of previous events, \code{N}, are carried
     forward assuming that no new events occur beyond censoring.
     If \code{last.value=FALSE}, then these variables are coded \code{NA}
     for easy identification allowing replacement with the desired values.
   }
   
   %% \item{probs}{ \code{tx.test} is generated from \code{x.train} in a
   %%   similar fashion to \code{tx.train}.  However, observations beyond
   %%   censoring contain \code{NA} for columns two and three containing
   %%   \eqn{v(t)} and \eqn{N(t-)} respectively.  In order to preform
   %%   prediction with \code{tx.test}, then these missing values need to be
   %%   handled accordingly.  We provide an automated method for this which
   %%   is returned as \code{X.base}.  It would be difficult to provide an
   %%   efficient method that randomly generates \eqn{v(t)} and \eqn{N(t-)}
   %%   since they are structurally dependent.  Rather, we
   %%   deterministically generate them in such a fashion that they define
   %%   a "middle" baseline, i.e. the median \eqn{N(t-)} as a function of
   %%   \eqn{t} for non-censored subjects (calculated from \code{tx.train});
   %%   the median is adjusted upward as the censoring increases with time
   %%   (if it does at all) to a maximum value controlled by the second
   %%   argument contained in \code{probs} since otherwise the "low" risk
   %%   subjects would be over-represented.  Each subject's previous value
   %%   is carried forward unless their computed sojourn time (which is a
   %%   function of the generated \eqn{N(t-)}) is longer than the the
   %%   median \eqn{v(t)} as a function of \eqn{t} (which is also
   %%   calculated from \code{tx.train}); the median is adjusted downward as
   %%   the censoring increases with time (if it does at all) to a maximum
   %%   value controlled by the first argument contained in \code{probs}. }

   %% \item{baseline}{experimental feature: still needs work}
}

\value{
   \code{recur.pre.bart} returns a list.
Besides the items listed below, the list has
 a \code{times} component giving the unique times and \code{K} which is the number of
 unique times.
   \item{y.train}{A vector of binary responses.}
   \item{tx.train}{A matrix with the rows of the training data.}
   \item{tx.test}{Generated from \code{x.train} (see discussion above included in
     the argument \code{last.value}).}
   %% \item{X.base}{Generated from \code{x.train} (see discussion above included in
   %%   the argument \code{probs}).}
   %% \item{sojourn}{Generated from \code{x.train} (see discussion above included in
   %%   the argument \code{probs}).}
   %% \item{pattern}{Generated from \code{x.train} (see discussion above included in
   %%   the argument \code{probs}).}
}
%% \references{
%% Sparapani, Rein, Tarima, Jackson, Meurer (2020).  Non-parametric recurrent
%% events analysis with BART and an application to the hospital admissions
%% of patients with diabetes. \emph{Biostatistics} doi:10.1093/biostatistics/kxy032
%% %%MCW Biostatistics Technical Report 64
%% %%\url{http://www.mcw.edu/FileLibrary/Groups/Biostatistics/TechReports/TechReports5175/tr064.pdf}

%% %% Xu, Daniels and Winterstein (2016).  Sequential BART for imputation of
%% %% missing covariates.  \emph{Biostatistics} doi:10.1093/biostatistics/kxw009
%% %% %% with supplement
%% %% \url{http://biostatistics.oxfordjournals.org/content/early/2016/03/15/biostatistics.kxw009/suppl/DC1}
%% }
%% \author{
%% Rodney Sparapani: \email{rsparapa@mcw.edu}
%% }
\seealso{
\code{\link{recur.bart}}
}
\examples{

data(bladder)
subset <- -which(bladder1$stop==0)
bladder0 <- bladder1[subset, ]
id <- unique(sort(bladder0$id))
N <- length(id)
L <- max(bladder0$enum)

times <- matrix(0, nrow=N, ncol=L)
dimnames(times)[[1]] <- paste0(id)

delta <- matrix(0, nrow=N, ncol=L)
dimnames(delta)[[1]] <- paste0(id)

x.train <- matrix(NA, nrow=N, ncol=3+2*L) ## add time-dependent cols too
dimnames(x.train)[[1]] <- paste0(id)
dimnames(x.train)[[2]] <- c('Pl', 'B6', 'Th', rep(c('number', 'size'), L))

for(i in 1:N) {
    h <- id[i]

    for(j in 1:L) {
        k <- which(bladder0$id==h & bladder0$enum==j)

        if(length(k)==1) {
            times[i, j] <- bladder0$stop[k]
            delta[i, j] <- (bladder0$status[k]==1)*1

            if(j==1) {
                x.train[i, 1] <- as.numeric(bladder0$treatment[k])==1
                x.train[i, 2] <- as.numeric(bladder0$treatment[k])==2
                x.train[i, 3] <- as.numeric(bladder0$treatment[k])==3
                x.train[i, 4] <- bladder0$number[k]
                x.train[i, 5] <- bladder0$size[k]
            }
            else if(delta[i, j]==1) {
                if(bladder0$rtumor[k]!='.')
                    x.train[i, 2*j+2] <- as.numeric(bladder0$rtumor[k])
                if(bladder0$rsize[k]!='.')
                    x.train[i, 2*j+3] <- as.numeric(bladder0$rsize[k])
            }
        }
    }
}

pre <- recur.pre.bart(times=times, delta=delta, x.train=x.train)

J <- nrow(pre$tx.train)
for(j in 1:J) {
    if(pre$tx.train[j, 3]>0) {
        pre$tx.train[j, 7] <- pre$tx.train[j, 7+pre$tx.train[j, 3]*2]
        pre$tx.train[j, 8] <- pre$tx.train[j, 8+pre$tx.train[j, 3]*2]
    }
}
pre$tx.train <- pre$tx.train[ , 1:8]

K <- pre$K
NK <- N*K
for(j in 1:NK) {
    if(pre$tx.test[j, 3]>0) {
        pre$tx.test[j, 7] <- pre$tx.test[j, 7+pre$tx.test[j, 3]*2]
        pre$tx.test[j, 8] <- pre$tx.test[j, 8+pre$tx.test[j, 3]*2]
    }
}
pre$tx.test <- pre$tx.test[ , 1:8]


## in bladder1 both number and size are recorded as integers
## from 1 to 8 however they are often missing for recurrences
## at baseline there are no missing and 1 is the mode of both
pre$tx.train[which(is.na(pre$tx.train[ , 7])), 7] <- 1
pre$tx.train[which(is.na(pre$tx.train[ , 8])), 8] <- 1
pre$tx.test[which(is.na(pre$tx.test[ , 7])), 7] <- 1
pre$tx.test[which(is.na(pre$tx.test[ , 8])), 8] <- 1

## it is a good idea to explore more sophisticated methods
## such as imputing the missing data with Sequential BART
## Xu, Daniels and Winterstein.  Sequential BART for imputation of missing
## covariates.  Biostatistics 2016 doi: 10.1093/biostatistics/kxw009
## http://biostatistics.oxfordjournals.org/content/early/2016/03/15/biostatistics.kxw009/suppl/DC1
## https://cran.r-project.org/package=sbart
## library(sbart)
## set.seed(21)
## train <- seqBART(xx=pre$tx.train, yy=NULL, datatype=rep(0, 6),
##                type=0, numskip=20, burn=1000)
## coarsen the imputed data same way as observed example data
## train$imputed5[which(train$imputed5[ , 7]<1), 7] <- 1
## train$imputed5[which(train$imputed5[ , 7]>8), 7] <- 8
## train$imputed5[ , 7] <- round(train$imputed5[ , 7])
## train$imputed5[which(train$imputed5[ , 8]<1), 8] <- 1
## train$imputed5[which(train$imputed5[ , 8]>8), 8] <- 8
## train$imputed5[ , 8] <- round(train$imputed5[ , 8])

## for Friedman's partial dependence, we need to estimate the whole cohort
## at each treatment assignment (and, average over those)
pre$tx.test <- rbind(pre$tx.test, pre$tx.test, pre$tx.test)
pre$tx.test[ , 4] <- c(rep(1, NK), rep(0, 2*NK))          ## Pl
pre$tx.test[ , 5] <- c(rep(0, NK), rep(1, NK), rep(0, NK))## B6
pre$tx.test[ , 6] <- c(rep(0, 2*NK), rep(1, NK))          ## Th

\dontrun{
## set.seed(99)
## post <- recur.bart(y.train=pre$y.train, x.train=pre$tx.train, x.test=pre$tx.test)
## depending on your performance, you may want to run in parallel if available
post <- mc.recur.bart(y.train=pre$y.train, x.train=pre$tx.train,
                      x.test=pre$tx.test, mc.cores=8, seed=99)

M <- nrow(post$yhat.test)
RI.B6.Pl <- matrix(0, nrow=M, ncol=K)
RI.Th.Pl <- matrix(0, nrow=M, ncol=K)
RI.Th.B6 <- matrix(0, nrow=M, ncol=K)

for(j in 1:K) {
    h <- seq(j, NK, K)
    RI.B6.Pl[ , j] <- apply(post$prob.test[ , h+NK]/
                            post$prob.test[ , h], 1, mean)
    RI.Th.Pl[ , j] <- apply(post$prob.test[ , h+2*NK]/
                            post$prob.test[ , h], 1, mean)
    RI.Th.B6[ , j] <- apply(post$prob.test[ , h+2*NK]/
                            post$prob.test[ , h+NK], 1, mean)
}

RI.B6.Pl.mu <- apply(RI.B6.Pl, 2, mean)
RI.B6.Pl.025 <- apply(RI.B6.Pl, 2, quantile, probs=0.025)
RI.B6.Pl.975 <- apply(RI.B6.Pl, 2, quantile, probs=0.975)

RI.Th.Pl.mu <- apply(RI.Th.Pl, 2, mean)
RI.Th.Pl.025 <- apply(RI.Th.Pl, 2, quantile, probs=0.025)
RI.Th.Pl.975 <- apply(RI.Th.Pl, 2, quantile, probs=0.975)

RI.Th.B6.mu <- apply(RI.Th.B6, 2, mean)
RI.Th.B6.025 <- apply(RI.Th.B6, 2, quantile, probs=0.025)
RI.Th.B6.975 <- apply(RI.Th.B6, 2, quantile, probs=0.975)

plot(post$times, RI.Th.Pl.mu, col='blue',
     log='y', main='Bladder cancer ex: Thiotepa vs. Placebo',
     type='l', ylim=c(0.1, 10), ylab='RI(t)', xlab='t (months)')
lines(post$times, RI.Th.Pl.025, col='red')
lines(post$times, RI.Th.Pl.975, col='red')
abline(h=1)

plot(post$times, RI.B6.Pl.mu, col='blue',
     log='y', main='Bladder cancer ex: Vitamin B6 vs. Placebo',
     type='l', ylim=c(0.1, 10), ylab='RI(t)', xlab='t (months)')
lines(post$times, RI.B6.Pl.025, col='red')
lines(post$times, RI.B6.Pl.975, col='red')
abline(h=1)

plot(post$times, RI.Th.B6.mu, col='blue',
     log='y', main='Bladder cancer ex: Thiotepa vs. Vitamin B6',
     type='l', ylim=c(0.1, 10), ylab='RI(t)', xlab='t (months)')
lines(post$times, RI.Th.B6.025, col='red')
lines(post$times, RI.Th.B6.975, col='red')
abline(h=1)

}
}
\keyword{data construction}
