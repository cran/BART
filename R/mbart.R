mbart = function(
x.train, y.train, x.test=matrix(0.0,0,0),
sparse=FALSE, theta=0, omega=1,
a=0.5, b=1, augment=FALSE, rho=NULL,
xinfo=matrix(0.0,0,0), usequants=FALSE,
cont=FALSE, rm.const=TRUE, ##tau.interval=0.95,
k=2.0, power=2.0, base=.95,
binaryOffset=NULL,
ntree=50L, numcut=100L,
ndpost=1000L, nskip=100L,
keepevery=1L,
nkeeptrain=ndpost, nkeeptest=ndpost,
#nkeeptestmean=ndpost,
nkeeptreedraws=ndpost,
printevery=100, transposed=FALSE
#treesaslists=FALSE
)
{

n = length(y.train)

## if(binaryOffset!=0)
##     stop('binaryOffset not supported by mbart')

if(length(binaryOffset)==0)
    binaryOffset <- qnorm(as.integer(table(y.train))/n)

if(!transposed) {
    temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                           cont=cont, xinfo=xinfo, rm.const=rm.const)
    x.train = t(temp$X)
    numcut = temp$numcut
    xinfo = temp$xinfo
    if(length(x.test)>0)
            x.test = t(bartModelMatrix(x.test[ , temp$rm.const]))
    rm.const <- temp$rm.const
    grp <- temp$grp
    rm(temp)
}
else {
    rm.const <- NULL
    grp <- NULL
}

if(n!=ncol(x.train))
    stop('The length of y.train and the number of rows in x.train must be identical')

p = nrow(x.train)
np = ncol(x.test)
if(length(rho)==0) rho <- p
if(length(rm.const)==0) rm.const <- 1:p
if(length(grp)==0) grp <- 1:p

## if(tau.interval>0.5) tau.interval=1-tau.interval

## tau=qlogis(1-0.5*tau.interval)/(k*sqrt(ntree))

#--------------------------------------------------
#set  nkeeps for thinning
if((nkeeptrain!=0) & ((ndpost %% nkeeptrain) != 0)) {
   nkeeptrain=ndpost
   cat('*****nkeeptrain set to ndpost\n')
}
if((nkeeptest!=0) & ((ndpost %% nkeeptest) != 0)) {
   nkeeptest=ndpost
   cat('*****nkeeptest set to ndpost\n')
}
## if((nkeeptestmean!=0) & ((ndpost %% nkeeptestmean) != 0)) {
##    nkeeptestmean=ndpost
##    cat('*****nkeeptestmean set to ndpost\n')
## }
if((nkeeptreedraws!=0) & ((ndpost %% nkeeptreedraws) != 0)) {
   nkeeptreedraws=ndpost
   cat('*****nkeeptreedraws set to ndpost\n')
}
C <- length(unique(sort(y.train)))
res = .Call("cmbart",
            n,  #number of observations in training data
            p,  #dimension of x
            np, #number of observations in test data
            x.train,   #p*n training data x
            as.integer(as.factor(y.train))-1, #0,1,,,,C-1 n*1 y
            C, #number of categories
            x.test,    #p*np test data x
            ntree,
            numcut,
            ndpost*keepevery,
            nskip,
            power,
            base,
            binaryOffset,
            3/(k*sqrt(ntree)),##tau,
            sparse,
            theta,
            omega,
            grp,
            a,
            b,
            rho,
            augment,
            nkeeptrain,
            nkeeptest,
#            nkeeptestmean,
            nkeeptreedraws,
            printevery,
            xinfo
)
res$C <- C

if(nkeeptrain>0) {
    res$prob.train <- pnorm(res$yhat.train)
    for(i in 1:n) {
        h <- (i-1)*C
        total <- apply(res$prob.train[ , h+1:C], 1, sum)
        for(j in 1:C) {
            res$prob.train[ , h+j] <- res$prob.train[ , h+j]/total
        }
    }
    res$prob.train.mean <- apply(res$prob.train, 2, mean)
    res$yhat.train.mean <- NULL
} else {
    res$yhat.train <- NULL
    res$yhat.train.mean <- NULL
}

if(np>0) {
    res$prob.test <- pnorm(res$yhat.test)
    for(i in 1:np) {
        h <- (i-1)*C
        total <- apply(res$prob.test[ , h+1:C], 1, sum)
        for(j in 1:C) {
            res$prob.test[ , h+j] <- res$prob.test[ , h+j]/total
        }
    }
    res$prob.test.mean <- apply(res$prob.test, 2, mean)
    res$yhat.test.mean <- NULL
} else {
    res$yhat.test <- NULL
    res$yhat.test.mean <- NULL
}

names. <- dimnames(x.train)[[1]]

if(nkeeptreedraws>0)
    names(res$treedraws$cutpoints) = names.

names. <- rep(names., each=C)

for(i in 1:p)
    for(j in 1:C) {
        h <- (i-1)*C+j
        names.[h] <- paste0(names.[h], '.', j)
    }

dimnames(res$varcount)[[2]] = as.list(names.)
dimnames(res$varprob)[[2]] = as.list(names.)
res$varcount.mean <- apply(res$varcount, 2, mean)
res$varprob.mean <- apply(res$varprob, 2, mean)
res$rm.const <- rm.const
res$binaryOffset=binaryOffset
attr(res, 'class')='mbart'
return(res)
}
