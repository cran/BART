library(BART)
library(MASS)

##options(figures='../vignettes/figures')

B <- getOption('mc.cores', 1)
figures = getOption('figures', default='NONE')

y = Boston$medv # median value
x.train = as.matrix(cbind(Boston[ , -c(5, 14)], Boston[ , 5]))
dimnames(x.train)[[2]][13] = 'nox'
N=length(y)   ## total sample size
post = mc.gbart(x.train, y, mc.cores=B, seed=99)

L=41
x=seq(min(x.train[ , 13]), max(x.train[ , 13]), length.out=L)

x.test = cbind(x.train[ , -13], x[1])
names(x.test)[13]='nox'
for(j in 2:L)
    x.test = rbind(x.test, cbind(x.train[ , -13], x[j]))

pred = predict(post, x.test, mc.cores=B)

partial = matrix(nrow=1000, ncol=L)
for(j in 1:L) {
    h=(j-1)*N+1:N
    partial[ , j] = apply(pred[ , h], 1, mean)
}

plot(x, apply(partial, 2, mean), type='l', lwd=2,
     ##xlab='nox', ylab='mdev',
     xlab='nox: Nitrogen Oxides air pollution',
     ylab='mdev: median home value (in thousands)',
     ylim=c(0, 50))
lines(x, apply(partial, 2, quantile, probs=0.025), lty=2, lwd=2)
lines(x, apply(partial, 2, quantile, probs=0.975), lty=2, lwd=2)
abline(h=c(0, 50), col='gray')
## uncomment for an ICE plot
## for(i in seq(1, N, by=10)) 
##     lines(x, apply(pred[ , seq(i, N*L, by=N)], 2, mean), type='l',
##           lty=3, col='green')

## model similar to that presented in Harrison and Rubinfeld (1978)
fit=lm(log(y)~I(rm^2)+age+I(log(dis))+I(log(rad))+tax+ptratio+
           I((black-0.63)^2)+I(log(lstat))+crim+zn+indus+chas+
           I((nox-0.55)^2), data=Boston)
summary(fit)
lines(x, mean(y)*exp((x-0.55)^2*fit$coefficients[14]), lty=3, col='red', lwd=2)

if(figures!='NONE')
    dev.copy2pdf(file=paste(figures, 'nox.pdf', sep='/'))

if(figures!='NONE')
    dev.copy2eps(file=paste(figures, 'nox.eps', sep='/'))
