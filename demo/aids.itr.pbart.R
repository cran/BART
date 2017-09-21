
library(BART)

data(ACTG175)

## exclude those who do not have CD4 count at 96 weeks
ex <- is.na(ACTG175$cd496)
table(ex)

## inclusion criteria are CD4 counts between 200 and 500
ACTG175$cd40 <- min(500, max(250, ACTG175$cd40))

## calculate relative CD4 decline
y <- ((ACTG175$cd496-ACTG175$cd40)/ACTG175$cd40)[!ex]
summary(y)

## 0=failure, 1=success
y <- 1*(y > -0.5)

## summarize CD4 outcomes
table(y, ACTG175$arms[!ex])

table(y, ACTG175$arms[!ex])/
    matrix(table(ACTG175$arms[!ex]), nrow=2, ncol=4, byrow=TRUE)

## drop unneeded and unwanted variables
## 1: 'pidnum' patient ID number
##14: 'str2' which will be handled by strat1 below
##15: 'strat' which will be handled by strat1-strat3 below
##17: 'treat' handled by arm0-arm3 below
##18: 'offtrt' indicator of off-treatment before 96 weeks
##20: 'cd420' CD4 T cell count at 20 weeks
##21: 'cd496' CD4 T cell count at 96 weeks
##22: 'r' missing CD4 T cell count at 96 weeks
##24: 'cd820' CD8 T cell count at 20 weeks
##25: 'cens' indicator of observing the event in days
##26: 'days' number of days until the primary endpoint
##27: 'arms' handled by arm0-arm3 below
train <- as.matrix(ACTG175)[!ex, -c(1, 14:15, 17, 18, 20:22, 24:27)]
train <- cbind(1*(ACTG175$strat[!ex]==1), 1*(ACTG175$strat[!ex]==2),
               1*(ACTG175$strat[!ex]==3), train)
dimnames(train)[[2]][1:3] <- paste0('strat', 1:3)
train <- cbind(1*(ACTG175$arms[!ex]==0), 1*(ACTG175$arms[!ex]==1),
               1*(ACTG175$arms[!ex]==2), 1*(ACTG175$arms[!ex]==3), train)
dimnames(train)[[2]][1:4] <- paste0('arm', 0:3)

N <- nrow(train)

test0 <- train; test0[ , 1:4] <- 0; test0[ , 1] <- 1
test1 <- train; test1[ , 1:4] <- 0; test1[ , 2] <- 1
test2 <- train; test2[ , 1:4] <- 0; test2[ , 3] <- 1
test3 <- train; test3[ , 1:4] <- 0; test3[ , 4] <- 1

test <- rbind(test0, test1, test2, test3)

set.seed(21)
post <- pbart(train, y, test)

## turn z-scores into probabilities
post$prob.test <- pnorm(post$yhat.test)

## average over the posterior samples
post$prob.test.mean <- apply(post$prob.test, 2, mean)

## place estimates for arms 0-3 next to each other for convenience
itr <- cbind(post$prob.test.mean[(1:N)], post$prob.test.mean[N+(1:N)],
             post$prob.test.mean[2*N+(1:N)], post$prob.test.mean[3*N+(1:N)])

## find the BART ITR for each patient
itr.pick <- integer(N)
for(i in 1:N) itr.pick[i] <- which(itr[i, ]==max(itr[i, ]))-1

## arms 0 and 3 (monotherapy) are never chosen
table(itr.pick)

## do arms 1 and 2 show treatment heterogeneity?
diff. <- apply(post$prob.test[ , 2*N+(1:N)]-post$prob.test[ , N+(1:N)], 2, mean)
plot(sort(diff.), type='h', main='ACTG175 trial: 50% CD4 decline from baseline at 96 weeks',
     xlab='Arm 2 (1) Preferable to the Right (Left)', ylab='Prob.Diff.: Arms 2 - 1')

library(rpart)
library(rpart.plot)

## make data frame for nicer names in the plot
var <- as.data.frame(train[ , -(1:4)])

dss <- rpart(diff. ~ var$age+var$gender+var$race+var$wtkg+var$cd40+var$cd80+
                 var$karnof+var$symptom+var$hemo+var$homo+var$drugs+var$z30+
                 var$zprior+var$oprior+var$strat1+var$strat2+var$strat3,
             method='anova', control=rpart.control(cp=0.1))
rpart.plot(dss, type=3, extra=101)

## if strat1==1 (antiretroviral naive), then arm 2 is better
## otherwise, arm 1
print(dss)

all0 <- apply(post$prob.test[ , (1:N)], 1, mean)
all1 <- apply(post$prob.test[ , N+(1:N)], 1, mean)
all2 <- apply(post$prob.test[ , 2*N+(1:N)], 1, mean)
all3 <- apply(post$prob.test[ , 3*N+(1:N)], 1, mean)

## BART ITR
BART.itr <- apply(post$prob.test[ , c(N+which(itr.pick==1), 2*N+which(itr.pick==2))], 1, mean)

test <- train
test[ , 1:4] <- 0
test[test[ , 5]==0, 2] <- 1
test[test[ , 5]==1, 3] <- 1

## BART ITR simple
BART.itr.simp <- pwbart(test, post$treedraws)
BART.itr.simp <- apply(pnorm(BART.itr.simp), 1, mean)

plot(density(BART.itr), xlab='Value', xlim=c(0.475, 0.775), lwd=2,
     main='ACTG175 trial: 50% CD4 decline from baseline at 96 weeks')
lines(density(BART.itr.simp), col='brown', lwd=2)
lines(density(all0), col='green', lwd=2)
lines(density(all1), col='red', lwd=2)
lines(density(all2), col='blue', lwd=2)
lines(density(all3), col='yellow', lwd=2)
legend('topleft', legend=c('All Arm 0 (ZDV only)', 'All Arm 1 (ZDV+DDI)',
                           'All Arm 2 (ZDV+DDC)', 'All Arm 3 (DDI only)',
                           'BART ITR simple', 'BART ITR'),
       col=c('green', 'red', 'blue', 'yellow', 'brown', 'black'), lty=1, lwd=2)
