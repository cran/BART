
library(BART)

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

##set.seed(99)
##post=mbart2(x.train, y.train, x.test, type='pbart')
post=mc.mbart2(x.train, y.train, x.test, type='pbart', mc.cores=8, seed=99)
## check=predict(post, x.test, mc.cores=8)
## print(cor(post$prob.test.mean, check$prob.test.mean)^2)

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
pred=predict(post, x.test, mc.cores=8)
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

k=1:K
size.test.LH1=double(2*K)
size.test.LH1[2*k-1]=size.test.025[k]
size.test.LH1[2*k]=size.test.975[k]
size.test.LH2=double(2*K)
size.test.LH2[2*k-1]=size.test.025[k+K]
size.test.LH2[2*k]=size.test.975[k+K]
plot(factor(k, labels=c('bird', 'fish', 'invert', 'other', 'reptile')),
     rep(1, K), col=k, type='n', lwd=2, lty=0,
             xlim=c(1, K), ylim=c(0, 0.6), ylab='Probability',
     sub="Multinomial BART\nFriedman's partial dependence function")
points(k-0.05, size.test.mean[k+K], lwd=2, col=1)
points(k+0.05, size.test.mean[k], lwd=2, col=2)
for(k in 1:K) {
    lines(rep(k+0.05, 2), size.test.LH1[c(2*k-1, 2*k)], lwd=2, col=2)
    lines(rep(k-0.05, 2), size.test.LH2[c(2*k-1, 2*k)], lwd=2, lty=2)
}
legend('topright', legend=c('Small', 'Large'),
        pch=1, col=1:2, lty=2:1, lwd=2)
dev.copy2pdf(file='../vignettes/figures/alligator.pdf')

## plot(factor(1:K, labels=c('bird', 'fish', 'invert', 'other', 'reptile')),
##      rep(1, K), col=1:K, type='n', lwd=2, lty=0,
##              xlim=c(1, K), ylim=c(0, 0.5), ylab='Prob.',
##      sub="Multinomial BART\nFriedman's partial dependence function")
## points(1:K, size.test.mean[1:K+K], lwd=2, col=1)
## lines(1:K, size.test.025[1:K+K], lwd=2, col=1, lty=2)
## lines(1:K, size.test.975[1:K+K], lwd=2, col=1, lty=2)
## points(1:K, size.test.mean[1:K], lwd=2, col=2)
## lines(1:K, size.test.025[1:K], lwd=2, col=2, lty=2)
## lines(1:K, size.test.975[1:K], lwd=2, col=2, lty=2)
## legend('topright', legend=c('Small', 'Large'),
##         pch=1, col=1:2)
##dev.copy2pdf(file='../vignettes/figures/alligator.pdf')
