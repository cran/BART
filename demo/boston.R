library(MASS)

B <- getOption('mc.cores', 1)
figures = getOption('figures', default='NONE')

x = Boston[,c(6,13)] #rm=number of rooms and lstat= percent lower status
y = Boston$medv # median value

head(cbind(x, y))

par(mfrow=c(2,2))
par(mai=c(.8,.8,.2,.2))
plot(x[,1],y,xlab="x1=rm",ylab="y=mdev",cex.axis=1.3,cex.lab=1.2)
plot(x[,2],y,xlab="x2=lstat",ylab="y=mdev",cex.axis=1.3,cex.lab=1.2)
plot(x[,1],x[,2],xlab="x1=rm",ylab="x2=lstat",cex.axis=1.3,cex.lab=1.2) > par(mfrow=c(1,1))
if(figures!='NONE')
    dev.copy2pdf(file=paste(figures, 'boston1.pdf', sep='/'))

library(BART) ## load library
set.seed(99)  ## MCMC posterior sampling: set seed for reproducibility
nd=200        ## number of draws to keep
burn=50       ## number of draws to discard
bf = wbart(x,y,nskip=burn,ndpost=nd)
plot(bf$sigma, ylab='post$sigma', type="l")
abline(v=burn,lwd=2,col="red")
if(figures!='NONE')
    dev.copy2pdf(file=paste(figures, 'boston2.pdf', sep='/'))

lmf = lm(y~.,data.frame(x,y))
fitmat = cbind(y,bf$yhat.train.mean,lmf$fitted.values)
colnames(fitmat)=c("y","BART","Linear")
cor(fitmat)
pairs(fitmat)
if(figures!='NONE')
    dev.copy2pdf(file=paste(figures, 'boston3.pdf', sep='/'))

ii = order(bf$yhat.train.mean) ## order observations by predicted value
boxplot(bf$yhat.train[,ii], ylab='post$yhat.train') ## boxplots of f(x) draws
if(figures!='NONE')
    dev.copy2pdf(file=paste(figures, 'boston4.pdf', sep='/'))

n=length(y)   ## total sample size
set.seed(14)  ## Dave Keon, greatest Leaf of all time!
## ii = sample(1:n,floor(.75*n)) ## indices for train data, 75% of data
## xtrain=x[ii,]; ytrain=y[ii]   ## training data
## xtest=x[-ii,]; ytest=y[-ii]   ## test data
## cat("train sample size is ",length(ytrain),"\n")
## cat("test sample size is ",length(ytest),"\n")
i <- sample(1:n, floor(0.75 * n))
xtrain <- x[i, ]; ytrain = y[i]
xtest <- x[-i, ]; ytest = y[-i]
cat("training sample size = ", length(ytrain), "\n")
cat("testing sample size = ", length(ytest), "\n")
set.seed(99)
bfp1 = wbart(xtrain,ytrain,xtest)
set.seed(99)
bfp2 = wbart(xtrain,ytrain)
yhat = predict(bfp2, as.matrix(xtest), mc.cores=B)
set.seed(4) #Bobby Orr's jersey number is the seed
bfthin = wbart(xtrain,ytrain,nskip=1000,ndpost=10000,
               nkeeptrain=0,nkeeptest=0,nkeeptestmean=0,nkeeptreedraws=200)
yhatthin = predict(bfthin, as.matrix(xtest), mc.cores=B)
fmat=cbind(ytest,bfp1$yhat.test.mean,apply(yhatthin,2,mean))
colnames(fmat) = c("y","yhat","yhatThin")
##colnames(fmat) = c("y","BARTpred","BARTpredThin")
pairs(fmat)
if(figures!='NONE')
    dev.copy2pdf(file=paste(figures, 'boston5.pdf', sep='/'))

y.train = ytrain
x.train = as.matrix(Boston[ii, -14])
set.seed(12) ## Aaron Rodgers
post4 = wbart(x.train, y.train)

N=379
L=41
x=seq(min(x.train[ , 13]), max(x.train[ , 13]), length.out=L)

x.test = cbind(x.train[ , -13], x[1])
names(x.test)[13]='lstat'
for(j in 2:L)
    x.test = rbind(x.test, cbind(x.train[ , -13], x[j]))

pred1 = predict(post4, x.test, mc.cores=B)

partial = matrix(nrow=1000, ncol=L)
for(j in 1:L) {
    h=(j-1)*N+1:N
    partial[ , j] = apply(pred1[ , h], 1, mean)
}

plot(x, apply(partial, 2, mean), type='l',
     xlab='lstat', ylab='mdev',
     ##xlab='lstat: percent lower status', ylab='mdev: median home value',
     ylim=c(10, 50))
lines(x, apply(partial, 2, quantile, probs=0.025), lty=2)
lines(x, apply(partial, 2, quantile, probs=0.975), lty=2)
if(figures!='NONE')
    dev.copy2pdf(file=paste(figures, 'boston6.pdf', sep='/'))

H = 6
x.train.mean = apply(x.train, 2, mean)
x.test = matrix(x.train.mean, nrow=H, ncol=13, byrow=TRUE)
dimnames(x.test)[[2]]=names(x.train.mean)
##lstat = range(x.train[ , 'lstat'])
lstat = quantile(x.train[ , 'lstat'], probs=c(0.05, 0.95))
lstat = seq(lstat[1], lstat[2], length.out=H)
lstat.diff = lstat[2]-lstat[1]
x.test[ , 'lstat'] = lstat

pred2 = predict(post4, x.test, mc.cores=B)

diff = (pred2[ , 2:H]-pred2[ , 1:(H-1)])/lstat.diff
diff.mean = apply(diff, 2, mean)
diff.025 = apply(diff, 2, quantile, probs=0.025)
diff.975 = apply(diff, 2, quantile, probs=0.975)

plot(lstat[-H], diff.mean, type='l',
     xlab='lstat', ylab='conditional effect',
     ylim=c(min(min(diff.025), -max(diff.975)),
            max(-min(diff.025), max(diff.975))))
lines(lstat[-H], diff.025, type='l', lty=2)
lines(lstat[-H], diff.975, type='l', lty=2)
abline(h=0, col='gray')

