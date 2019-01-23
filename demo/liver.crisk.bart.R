
library(BART)

B <- getOption('mc.cores', 1)
figures = getOption('figures', default='NONE')

data(transplant)

pfit <- survfit(Surv(futime, event) ~ abo, transplant)

# competing risks for type O
plot(pfit[4,], xscale=7, xmax=735, col=1:3, lwd=2, ylim=c(0, 1),
       xlab='t (weeks)', ylab='Aalen-Johansen (AJ) CI(t)')
    legend(450, .4, c("Death", "Transplant", "Withdrawal"), col=1:3, lwd=2)
## plot(pfit[4,], xscale=30.5, xmax=735, col=1:3, lwd=2, ylim=c(0, 1),
##        xlab='t (months)', ylab='Aalen-Johansen (AJ) CI(t)')
##     legend(450, .4, c("Death", "Transplant", "Withdrawal"), col=1:3, lwd=2)

delta <- (as.numeric(transplant$event)-1)
## recode so that delta=1 is cause of interest; delta=2 otherwise
delta[delta==1] <- 4
delta[delta==2] <- 1
delta[delta>1] <- 2
table(delta, transplant$event)

times <- pmax(1, ceiling(transplant$futime/7)) ## weeks
##times <- pmax(1, ceiling(transplant$futime/30.5)) ## months
table(times)

typeO <- 1*(transplant$abo=='O')
typeA <- 1*(transplant$abo=='A')
typeB <- 1*(transplant$abo=='B')
typeAB <- 1*(transplant$abo=='AB')
table(typeA, typeO)

x.train <- cbind(typeO, typeA, typeB, typeAB)

x.test <- cbind(1, 0, 0, 0)
dimnames(x.test)[[2]] <- dimnames(x.train)[[2]]

## run one long MCMC chain in one process
## set.seed(99)
## post <- crisk.bart(x.train=x.train, times=times, delta=delta, x.test=x.test)

## in the interest of time, consider speeding it up by parallel processing
## run "mc.cores" number of shorter MCMC chains in parallel processes
post <- mc.crisk.bart(x.train=x.train, times=times, delta=delta,
                      x.test=x.test, seed=99, mc.cores=B)

K <- post$K

typeO.cif.mean <- apply(post$cif.test, 2, mean)
typeO.cif.025 <- apply(post$cif.test, 2, quantile, probs=0.025)
typeO.cif.975 <- apply(post$cif.test, 2, quantile, probs=0.975)

plot(pfit[4,], xscale=7, xmax=735, col=1:3, lwd=2, ylim=c(0, 0.8),
       xlab='t (weeks)', ylab='CI(t)')
points(c(0, post$times)*7, c(0, typeO.cif.mean), col=4, type='s', lwd=2)
points(c(0, post$times)*7, c(0, typeO.cif.025), col=4, type='s', lwd=2, lty=2)
points(c(0, post$times)*7, c(0, typeO.cif.975), col=4, type='s', lwd=2, lty=2)
     legend(450, .4, c("Transplant(BART)", "Transplant(AJ)",
                       "Death(AJ)", "Withdrawal(AJ)"),
            col=c(4, 2, 1, 3), lwd=2)
if(figures!='NONE')
    dev.copy2pdf(file=paste(figures, 'liver-BART.pdf', sep='/'))

## plot(pfit[4,], xscale=30.5, xmax=735, col=1:3, lwd=2, ylim=c(0, 0.8),
##        xlab='t (months)', ylab='CI(t)')
## points(c(0, post$times)*30.5, c(0, typeO.cif.mean), col=4, type='s', lwd=2)
## points(c(0, post$times)*30.5, c(0, typeO.cif.025), col=4, type='s', lwd=2, lty=2)
## points(c(0, post$times)*30.5, c(0, typeO.cif.975), col=4, type='s', lwd=2, lty=2)
##      legend(450, .4, c("Transplant(BART)", "Transplant(AJ)",
##                        "Death(AJ)", "Withdrawal(AJ)"),
##             col=c(4, 2, 1, 3), lwd=2)

## check <- predict(post, newdata=post$tx.test, newdata2=post$tx.test2,
##                  mc.cores=B)

## print(c(post$surv.test.mean[1], check$surv.test.mean[1],
##         post$surv.test.mean[1]-check$surv.test.mean[1]), digits=22)

## print(all(round(post$surv.test.mean, digits=9)==
##     round(check$surv.test.mean, digits=9)))

## print(c(post$cif.test.mean[1], check$cif.test.mean[1],
##         post$cif.test.mean[1]-check$cif.test.mean[1]), digits=22)

## print(all(round(post$cif.test.mean, digits=9)==
##     round(check$cif.test.mean, digits=9)))

## print(c(post$cif.test2.mean[1], check$cif.test2.mean[1],
##         post$cif.test2.mean[1]-check$cif.test2.mean[1]), digits=22)

## print(all(round(post$cif.test2.mean, digits=9)==
##     round(check$cif.test2.mean, digits=9)))

## typeO.cif.mean <- apply(check$cif.test, 2, mean)
## typeO.cif.025 <- apply(check$cif.test, 2, quantile, probs=0.025)
## typeO.cif.975 <- apply(check$cif.test, 2, quantile, probs=0.975)

## plot(pfit[4,], xscale=7, xmax=735, col=1:3, lwd=2, ylim=c(0, 0.8),
##        xlab='t (weeks)', ylab='CI(t)')
## points(c(0, post$times)*7, c(0, typeO.cif.mean), col=4, type='s', lwd=2)
## points(c(0, post$times)*7, c(0, typeO.cif.025), col=4, type='s', lwd=2, lty=2)
## points(c(0, post$times)*7, c(0, typeO.cif.975), col=4, type='s', lwd=2, lty=2)
##      legend(450, .4, c("Transplant(BART)", "Transplant(AJ)",
##                        "Death(AJ)", "Withdrawal(AJ)"),
##             col=c(4, 2, 1, 3), lwd=2)

## cor(post$cif.test.mean, check$cif.test.mean)
## plot(post$cif.test.mean, check$cif.test.mean)
## abline(a=0, b=1)

## cor(post$cif.test2.mean, check$cif.test2.mean)
## plot(post$cif.test2.mean, check$cif.test2.mean)
## abline(a=0, b=1)
