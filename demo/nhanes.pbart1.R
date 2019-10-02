
library(BART)

B <- getOption('mc.cores', 1)
figures = getOption('figures', default='NONE')

data(arq)
str(arq)
arth <- as.matrix(arq)

N <- length(arth[ , 'riagendr'])
table(arth[ , 'riagendr'])
summary(arth[ , 'bmxbmi'])

post1 <- mc.pbart(x.train=arth[ , 5:10], y.train=arth[ , 4],
                 mc.cores=B, seed=99)

post2 <- mc.pbart(x.train=arth[ , 5:10], y.train=arth[ , 3],
                 mc.cores=B, seed=99)

bmxbmi <- seq(15, 45, by=5)
H <- length(bmxbmi)

for(i in 1:2)
    for(j in 1:H) {
        x. <- arth[ , 5:10]
        x.[ , 'riagendr'] <- i
        x.[ , 'bmxbmi'] <- bmxbmi[j]
        if(i==1 && j==1) x.test <- x.
        else x.test <- rbind(x.test, x.)
    }

table(x.test[ , 'riagendr'])
table(x.test[ , 'bmxbmi'])

pred1 <- predict(post1, newdata=x.test, mc.cores=B)
pred2 <- predict(post2, newdata=x.test, mc.cores=B)

M <- nrow(pred1$prob.test)
##Friedman's partial dependence function
pd1 <- matrix(nrow=M, ncol=H)
pd2 <- matrix(nrow=M, ncol=H)

par(mfrow=c(1, 2))

for(j in 1:H) {
    h <- (j-1)*N
    pd1[ , j] <- apply(pred1$prob.test[ , h+1:N], 1, mean)
    h <- h+N*H
    pd2[ , j] <- apply(pred1$prob.test[ , h+1:N], 1, mean)
}

pd1.mean <- apply(pd1, 2, mean)
pd2.mean <- apply(pd2, 2, mean)
pd1.025 <- apply(pd1, 2, quantile, probs=0.025)
pd2.025 <- apply(pd2, 2, quantile, probs=0.025)
pd1.975 <- apply(pd1, 2, quantile, probs=0.975)
pd2.975 <- apply(pd2, 2, quantile, probs=0.975)

plot(bmxbmi, pd1.mean, type='l', col='blue',
     ylim=0:1, xlab='BMI', ylab=expression(p(x)),
     sub='Low-back pain: M(blue) vs. F(red)')
     ##sub='Low-back/buttock pain: M(blue) vs. F(red)')
lines(bmxbmi, pd1.025, type='l', col='blue', lty=2)
lines(bmxbmi, pd1.975, type='l', col='blue', lty=2)
lines(bmxbmi, pd2.mean, type='l', col='red')
lines(bmxbmi, pd2.025, type='l', col='red', lty=2)
lines(bmxbmi, pd2.975, type='l', col='red', lty=2)
lines(bmxbmi, rep(0, H), type='l')
lines(bmxbmi, rep(1, H), type='l')

for(j in 1:H) {
    h <- (j-1)*N
    pd1[ , j] <- apply(pred2$prob.test[ , h+1:N], 1, mean)
    h <- h+N*H
    pd2[ , j] <- apply(pred2$prob.test[ , h+1:N], 1, mean)
}

pd1.mean <- apply(pd1, 2, mean)
pd2.mean <- apply(pd2, 2, mean)
pd1.025 <- apply(pd1, 2, quantile, probs=0.025)
pd2.025 <- apply(pd2, 2, quantile, probs=0.025)
pd1.975 <- apply(pd1, 2, quantile, probs=0.975)
pd2.975 <- apply(pd2, 2, quantile, probs=0.975)

plot(bmxbmi, pd1.mean, type='l', col='blue',
     ylim=0:1, xlab='BMI', ylab=expression(p(x)),
     sub='Neck pain: M(blue) vs. F(red)')
lines(bmxbmi, pd1.025, type='l', col='blue', lty=2)
lines(bmxbmi, pd1.975, type='l', col='blue', lty=2)
lines(bmxbmi, pd2.mean, type='l', col='red')
lines(bmxbmi, pd2.025, type='l', col='red', lty=2)
lines(bmxbmi, pd2.975, type='l', col='red', lty=2)
lines(bmxbmi, rep(0, H), type='l')
lines(bmxbmi, rep(1, H), type='l')

par(mfrow=c(1, 1))

if(figures!='NONE')
    dev.copy2pdf(file=paste(figures, 'chronic-pain1.pdf', sep='/'))
##dev.copy2pdf(file='../vignettes/figures/chronic-pain1.pdf')

