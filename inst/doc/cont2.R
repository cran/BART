## ----setup, echo=FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----colfun,echo=FALSE---------------------------------------------------
library(knitr)
#Color Format
colFmt = function(x,color){
##   outputFormat = opts_knit$get("rmarkdown.pandoc.to")
     paste("\\textcolor{",color,"}{",x,"}",sep="")
}

## ----s1-1, include=TRUE, echo=TRUE,cache=TRUE----------------------------
library(MASS)
x = Boston[,c(6,13)] #rm=number of rooms and lstat= percent lower status
y = Boston$medv # median value
head(cbind(x,y))

## ----pl-dat, include=TRUE, echo=TRUE,out.width='80%',fig.align='center', dependson="s1-1"----
par(mfrow=c(2,2))
par(mai=c(.8,.8,.2,.2))
plot(x[,1],y,xlab="x1=rm",ylab="y=mdev",cex.axis=1.3,cex.lab=1.2)
plot(x[,2],y,xlab="x2=lstat",ylab="y=mdev",cex.axis=1.3,cex.lab=1.2)
plot(x[,1],x[,2],xlab="x1=rm",ylab="x2=lstat",cex.axis=1.3,cex.lab=1.2)

## ----s1-2, include=TRUE, echo=TRUE,cache=TRUE,dependson="s1-1",message=FALSE----
library(BART) #BART package
set.seed(99) #MCMC, so set the seed
nd=200 # number of kept draws
burn=50 # number of burn in draws
bf = wbart(x,y,nskip=burn,ndpost=nd)

## ----example, include=TRUE, echo=TRUE,dependson="s1-2",collapse=TRUE-----
names(bf)
length(bf$sigma)
length(bf$yhat.train.mean)
dim(bf$yhat.train)

## ----plotsigma, include=TRUE, echo=TRUE,dependson="s1-2", out.width='50%',fig.align='center'----
plot(bf$sigma)
abline(v=burn,lwd=2,col="red")

## ----comparison, include=TRUE, echo=TRUE,dependson="s1-2", out.width='50%',fig.align='center'----
lmf = lm(y~.,data.frame(x,y))
fitmat = cbind(y,bf$yhat.train.mean,lmf$fitted.values)
colnames(fitmat)=c("y","BART","Linear")
cor(fitmat)
pairs(fitmat)

## ----boxplots, include=TRUE, echo=TRUE,out.width='70%',fig.align='center',dependson="s1-2"----
ii = order(bf$yhat.train.mean) #order observations by predicted value
boxplot(bf$yhat.train[,ii]) #boxplots of f(x) draws

## ----ttsplit, include=TRUE, echo=TRUE,cache=TRUE,dependson="s1-1"--------
n=length(y) #total sample size
set.seed(14)  # Dave Keon, greatest Leaf of all time!
ii = sample(1:n,floor(.75*n)) # indices for train data, 75% of data
xtrain=x[ii,]; ytrain=y[ii] # training data
xtest=x[-ii,]; ytest=y[-ii] # test data
cat("train sample size is ",length(ytrain)," and test sample size is ",length(ytest),"\n")

## ----pred1,include=FALSE, echo=TRUE,cache=TRUE,message=FALSE,dependson="ttsplit"----
set.seed(99)
bfp1 = wbart(xtrain,ytrain,xtest) #predict.wbart wants a matrix

## ----output1, include=TRUE, echo=TRUE,dependson="pred1",collapse=TRUE----
dim(bfp1$yhat.test)
length(bfp1$yhat.test.mean)

## ----pred2,include=TRUE, echo=TRUE,cache=TRUE,message=FALSE,dependson="ttsplit"----
set.seed(99)
bfp2 = wbart(xtrain,ytrain)
yhat = predict(bfp2,as.matrix(xtest)) #predict wants a matrix

## ----output2, include=TRUE, echo=TRUE,dependson="pred2",collapse=TRUE----
dim(yhat)
summary(as.double(yhat-bfp1$yhat.test))

## ----bfthin,include=TRUE, echo=TRUE,cache=TRUE,message=FALSE,dependson="ttsplit"----
set.seed(4) #Bobby Orr's jersey number is the seed
bfthin = wbart(xtrain,ytrain,nskip=1000,ndpost=10000,
                     nkeeptrain=0,nkeeptest=0,nkeeptestmean=0,nkeeptreedraws=200)
yhatthin = predict(bfthin,as.matrix(xtest)) #predict wants a matrix

## ----output3, include=TRUE, echo=TRUE,collapse=TRUE----------------------
dim(bfthin$yhat.train)
dim(yhatthin)

## ----output4, include=TRUE, echo=TRUE,out.width='60%',fig.align='center',dependson="ttsplit"----
fmat=cbind(ytest,bfp1$yhat.test.mean,apply(yhatthin,2,mean))
colnames(fmat) = c("y","BARTpred","BARTpredThin")
pairs(fmat)

