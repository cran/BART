\name{alligator}
\alias{alligator}
\docType{data}
\title{ American alligator Food Choice }

\description{ In 1985, American alligators were harvested by hunters
from August 26 to September 30 in peninsular Florida from lakes Oklawaha
(Putnam County), George (Putnam and Volusia counties), Hancock (Polk
County) and Trafford (Collier County). Lake, length and sex were
recorded for each alligator. Stomachs from a sample of alligators
1.09-3.89m long were frozen prior to analysis. After thawing, stomach
contents were removed and separated and food items were identified and
tallied. Volumes were determined by water displacement. The stomach
contents of 219 alligators were classified into five categories of
primary food choice: Fish
(the most common primary food choice), Invertebrate (snails, insects,
crayfish, etc.), Reptile (turtles, alligators), Bird, and Other
(amphibians, plants, household pets, stones, and other debris).  }

\usage{data(alligator)}
\format{
  A data frame with 80 observations on the following 5 variables.
  \describe{
    \item{\code{lake}}{a factor with levels \code{George} \code{Hancock} \code{Oklawaha} \code{Trafford}}
    \item{\code{sex}}{a factor with levels \code{female} \code{male}}
    \item{\code{size}}{alligator size, a factor with levels \code{large} (>2.3m) \code{small} (<=2.3m)}
    \item{\code{food}}{primary food choice, a factor with levels \code{bird} \code{fish} \code{invert} \code{other} \code{reptile}}
    \item{\code{count}}{cell frequency, a numeric vector}
  }
}

\details{ The table contains a fair number of 0 counts.  \code{food} is
the response variable.  \code{fish} is the most frequent choice, and
often taken as a baseline category in multinomial response models.  }

\source{ 
  Agresti, A. (2002).
  \emph{Categorical Data Analysis},
  New York: Wiley, 2nd Ed., Table 7.1
}

\references{
  Delany MF, Linda SB, Moore CT (1999).
  "Diet and condition of American alligators in 4 Florida lakes."
  In \emph{Proceedings of the Annual Conference of the Southeastern
    Association of Fish and Wildlife Agencies}, \bold{53},
  375--389. 
}

\examples{

data(alligator)
     
## nnet::multinom Multinomial logit model fit with neural nets
fit <- multinom(food ~ lake+size+sex, data=alligator, weights=count)

summary(fit$fitted.values)
## 1=bird, 2=fish, 3=invert, 4=other, 5=reptile

(L=length(alligator$count))
(N=sum(alligator$count))
y.train=integer(N)
x.train=matrix(nrow=N, ncol=3)
x.test=matrix(nrow=L, ncol=3)
k=1
for(i in 1:L) {
    x.test[i, ]=as.integer(
        c(alligator$lake[i], alligator$size[i], alligator$sex[i]))
    if(alligator$count[i]>0)
        for(j in 1:alligator$count[i]) {
            y.train[k]=as.integer(alligator$food[i])
            x.train[k, ]=as.integer(
                c(alligator$lake[i], alligator$size[i], alligator$sex[i]))
            k=k+1
        }
}
table(y.train)
##test mbart with token run to ensure installation works
set.seed(99)
check = mbart(x.train, y.train, nskip=1, ndpost=1)

\dontrun{
set.seed(99)
check = mbart(x.train, y.train, nskip=1, ndpost=1)
post=mbart(x.train, y.train, x.test)

##post=mc.mbart(x.train, y.train, x.test, mc.cores=8, seed=99)
##check=predict(post, x.test, mc.cores=8)
##print(cor(post$prob.test.mean, check$prob.test.mean)^2)

par(mfrow=c(3, 2))
K=5
for(j in 1:5) {
    h=seq(j, L*K, K)
    print(cor(fit$fitted.values[ , j], post$prob.test.mean[h])^2)
    plot(fit$fitted.values[ , j], post$prob.test.mean[h],
         xlim=0:1, ylim=0:1,
         xlab=paste0('NN: Est. Prob. j=', j),
         ylab=paste0('BART: Est. Prob. j=', j))
    abline(a=0, b=1)
}
par(mfrow=c(1, 1))

L=16
x.test=matrix(nrow=L, ncol=3)
k=1
for(size in 1:2)
    for(sex in 1:2)
        for(lake in 1:4) {
            x.test[k, ]=c(lake, size, sex)
            k=k+1
        }
x.test

## two sizes: 1=large: >2.3m, 2=small: <=2.3m
pred=predict(post, x.test)
##pred=predict(post, x.test, mc.cores=8)
ndpost=nrow(pred$prob.test)

size.test=matrix(nrow=ndpost, ncol=K*2)
for(i in 1:K) {
    j=seq(i, L*K/2, K) ## large
    size.test[ , i]=apply(pred$prob.test[ , j], 1, mean)
    j=j+L*K/2 ## small
    size.test[ , i+K]=apply(pred$prob.test[ , j], 1, mean)
}
size.test.mean=apply(size.test, 2, mean)
size.test.025=apply(size.test, 2, quantile, probs=0.025)
size.test.975=apply(size.test, 2, quantile, probs=0.975)

plot(factor(1:K, labels=c('bird', 'fish', 'invert', 'other', 'reptile')),
     rep(1, K), col=1:K, type='n', lwd=1, lty=0,
             xlim=c(1, K), ylim=c(0, 0.5), ylab='Prob.',
     sub="Multinomial BART\nFriedman's partial dependence function")
points(1:K, size.test.mean[1:K+K], col=1)
lines(1:K, size.test.025[1:K+K], col=1, lty=2)
lines(1:K, size.test.975[1:K+K], col=1, lty=2)
points(1:K, size.test.mean[1:K], col=2)
lines(1:K, size.test.025[1:K], col=2, lty=2)
lines(1:K, size.test.975[1:K], col=2, lty=2)
## legend('topright', legend=c('Small', 'Large'),
##        pch=1, col=1:2)

}
}
\keyword{datasets}
