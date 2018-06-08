
library(BART)

data(arq)
str(arq)
arth <- as.matrix(arq)

N <- length(arth[ , 'riagendr'])
table(arth[ , 'riagendr'])
summary(arth[ , 'bmxbmi'])

post <- mc.pbart(x.train=arth[ , 5:10], y.train=arth[ , 4],
                 mc.cores=8, seed=99)

bmxbmi <- seq(15, 85, by=5)

for(i in 1:2) for(j in 1:15) {
                  x. <- arth[ , 5:10]
                  x.[ , 'riagendr'] <- i
                  x.[ , 'bmxbmi'] <- bmxbmi[j]
                  if(i==1 && j==1) x.test <- x.
                  else x.test <- rbind(x.test, x.)
              }

table(x.test[ , 'riagendr'])
table(x.test[ , 'bmxbmi'])

pred <- predict(post, newdata=x.test, mc.cores=8)

M <- nrow(pred$prob.test)

##Friedman's partial dependence function
pd1 <- matrix(nrow=M, ncol=15)
pd2 <- matrix(nrow=M, ncol=15)
for(j in 1:15) {
    h <- (j-1)*N
    pd1[ , j] <- apply(pred$prob.test[ , h+1:N], 1, mean)
    h <- h+N*15
    pd2[ , j] <- apply(pred$prob.test[ , h+1:N], 1, mean)
}

pd1.mean <- apply(pd1, 2, mean)
pd2.mean <- apply(pd2, 2, mean)
pd1.025 <- apply(pd1, 2, quantile, probs=0.025)
pd2.025 <- apply(pd2, 2, quantile, probs=0.025)
pd1.975 <- apply(pd1, 2, quantile, probs=0.975)
pd2.975 <- apply(pd2, 2, quantile, probs=0.975)

plot(bmxbmi, pd1.mean, type='l', col='blue',
     ylim=0:1, xlab='BMI', ylab=expression(Phi(f(x))),
     sub='Unweighted NHANES chronic low-back/buttock pain: M(blue) vs. F(red)')
## lines(bmxbmi, pd1.025, type='l', col='blue', lty=2)
## lines(bmxbmi, pd1.975, type='l', col='blue', lty=2)
lines(bmxbmi, pd2.mean, type='l', col='red')
## lines(bmxbmi, pd2.025, type='l', col='red', lty=2)
## lines(bmxbmi, pd2.975, type='l', col='red', lty=2)
##dev.copy2pdf(file='../vignettes/figures/clbp.pdf')

##incorporate survey weights into the posterior
wt.pd1 <- matrix(nrow=M, ncol=15)
wt.pd2 <- matrix(nrow=M, ncol=15)
for(j in 1:15) {
    h <- (j-1)*N
    wt.pd1[ , j] <- pred$prob.test[ , h+1:N] %*% (arth[ , 'wtint2yr']/sum(arth[ , 'wtint2yr']))
    h <- h+N*15
    wt.pd2[ , j] <- pred$prob.test[ , h+1:N] %*% (arth[ , 'wtint2yr']/sum(arth[ , 'wtint2yr']))
}

wt.pd1.mean <- apply(wt.pd1, 2, mean)
wt.pd2.mean <- apply(wt.pd2, 2, mean)
plot(bmxbmi, wt.pd1.mean, type='l', col='blue',
     ylim=0:1, xlab='BMI', ylab=expression(Phi(f(x))),
     sub='Weighted NHANES chronic low-back/buttock pain: M(blue) vs. F(red)')
lines(bmxbmi, wt.pd2.mean, type='l', col='red')

plot(bmxbmi, pd1.mean, type='l', col='blue',
     ylim=c(0.12, 0.36), xlab='BMI', ylab=expression(Phi(f(x))),
     sub='NHANES chronic low-back/buttock pain: M(blue) vs. F(red)')
lines(bmxbmi, pd2.mean, type='l', col='red')
lines(bmxbmi, wt.pd1.mean, type='l', col='blue', lty=2)
lines(bmxbmi, wt.pd2.mean, type='l', col='red', lty=2)
