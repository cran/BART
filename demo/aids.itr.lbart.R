
library(BART)

data(ACTG175)

## exclude those who do not have CD4 count at 96 weeks
## or their baseline CD4 is ineligible
## (inclusion criteria: CD4 counts between 200 and 500)
ex <- is.na(ACTG175$cd496) | ACTG175$cd40<200 | ACTG175$cd40>500
table(ex)

## calculate relative CD4 decline
y <- ((ACTG175$cd496-ACTG175$cd40)/ACTG175$cd40)[!ex]
summary(y)

## 0=failure, 1=success
y <- 1*(y > -0.5)

## summarize CD4 outcomes
a <- table(y, ACTG175$arms[!ex])
b <- apply(a, 2, sum)
a/rbind(b, b)

## drop unneeded and unwanted variables
## 1: 'pidnum' patient ID number
##10: 'zprior' zidovudine use prior to treatment initiation
##14: 'str2' which will be handled by strat1 below
##15: 'strat' which will be handled by strat1-strat3 below
##17: 'treat' handled by arm0-arm3 below
##18: 'offtrt' indicator of off-treatment before 96 weeks
##20: 'cd420' CD4 T cell count at 20 weeks
##21: 'cd496' CD4 T cell count at 96 weeks
##22: 'r' missing CD4 T cell count at 96 weeks
##23: 'cd80' CD8 T cell count at baseline 
##24: 'cd820' CD8 T cell count at 20 weeks
##25: 'cens' indicator of observing the event in days
##26: 'days' number of days until the primary endpoint
##27: 'arms' handled by arm0-arm3 below
train <- as.matrix(ACTG175)[!ex, -c(1, 10, 14:15, 17, 18, 20:27)]
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
post <- mc.lbart(train, y, test, mc.cores=8)

## place estimates for arms 0-3 next to each other for convenience
itr <- cbind(post$prob.test.mean[(1:N)],
             post$prob.test.mean[N+(1:N)],
             post$prob.test.mean[2*N+(1:N)],
             post$prob.test.mean[3*N+(1:N)])

## find the BART ITR for each patient
itr.pick <- integer(N)
for(i in 1:N) itr.pick[i] <- which(itr[i, ]==max(itr[i, ]))-1

## arms 0 and 3 (monotherapy) are never chosen
table(itr.pick)

## do arms 1 and 2 show treatment heterogeneity?
diff. <- apply(post$prob.test[ , 2*N+(1:N)]-
               post$prob.test[ , N+(1:N)], 2, mean)

plot(sort(diff.), type='h',
     main='ACTG175 trial: 50% CD4 decline from baseline at 96 weeks',
     xlab='Arm 2 (1) Preferable to the Right (Left)',
     ylab='Prob.Diff.: Arms 2 - 1')

library(rpart)
library(rpart.plot)

## make data frame for nicer names in the plot
var <- as.data.frame(train[ , -(1:4)])
##check <- complete.cases(var)
##table(check)

dss <- rpart(diff. ~ var$age+var$gender+var$race+var$wtkg+##var$cd80+
                 var$karnof+var$symptom+var$hemo+var$homo+var$drugs+var$z30+
                 var$oprior+var$strat1+var$strat2+var$strat3,
             method='anova', control=rpart.control(cp=0.05))
rpart.plot(dss, type=3, extra=101)

## if drugs==1, then arm 1
## otherwise, arm 2
print(dss)

all0 <- apply(post$prob.test[ , (1:N)], 1, mean)
all1 <- apply(post$prob.test[ , N+(1:N)], 1, mean)
all2 <- apply(post$prob.test[ , 2*N+(1:N)], 1, mean)
all3 <- apply(post$prob.test[ , 3*N+(1:N)], 1, mean)

## BART ITR
BART.itr <- apply(post$prob.test[ , c(N+which(itr.pick==1), 2*N+which(itr.pick==2))], 1, mean)

test <- train
test[ , 1:4] <- 0
test[var$drugs==1, 2] <- 1
test[ , 3] <- 1-test[ , 2]

table(test[ , 2], test[ , 3])

## BART ITR simple
BART.itr.simp <- predict(post, newdata=test, mc.cores=8)
BART.itr.simp$prob <- apply(BART.itr.simp$prob.test, 1, mean)

plot(density(BART.itr), xlab='Value', lwd=2, xlim=c(0.7, 0.95), 
     main='ACTG175 trial: 50% CD4 decline from baseline at 96 weeks')
lines(density(BART.itr.simp$prob), col='brown', lwd=2)
lines(density(all0), col='green', lwd=2)
lines(density(all1), col='red', lwd=2)
lines(density(all2), col='blue', lwd=2)
lines(density(all3), col='yellow', lwd=2)
legend('topleft', legend=c('All Arm 0 (ZDV only)',
                           'All Arm 1 (ZDV+DDI)',
                           'All Arm 2 (ZDV+DDC)',
                           'All Arm 3 (DDI only)',
                           'BART ITR simple',
                           'BART ITR'),
       col=c('green', 'red', 'blue', 'yellow', 'brown', 'black'),
       lty=1, lwd=2)
##dev.copy2pdf(file='aids.itr.lbart.pdf')
