
library(BART)

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

times <- pmax(1, transplant$futime/7) ## weeks
##times <- pmax(1, ceiling(transplant$futime/30.5)) ## months

typeO <- 1*(transplant$abo=='O')
typeA <- 1*(transplant$abo=='A')
typeB <- 1*(transplant$abo=='B')
typeAB <- 1*(transplant$abo=='AB')
table(typeA, typeO)

x.train <- cbind(typeO, typeA, typeB, typeAB)

x.test <- cbind(1, 0, 0, 0)
dimnames(x.test)[[2]] <- dimnames(x.train)[[2]]

## pre <- crisk.pre.bart(x.train=x.train, times=times, delta=delta,
##                       x.test=x.test, K=50)

## run one long MCMC chain in one process
## set.seed(99)
## post <- crisk.bart(x.train=x.train, times=times, delta=delta, x.test=x.test)

## in the interest of time, consider speeding it up by parallel processing
## run "mc.cores" number of shorter MCMC chains in parallel processes
post <- mc.crisk.bart(x.train=x.train, times=times, delta=delta, x.test=x.test,
                      K=50, seed=99, mc.cores=8)

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
## plot(pfit[4,], xscale=30.5, xmax=735, col=1:3, lwd=2, ylim=c(0, 0.8),
##        xlab='t (months)', ylab='CI(t)')
## points(c(0, post$times)*30.5, c(0, typeO.cif.mean), col=4, type='s', lwd=2)
## points(c(0, post$times)*30.5, c(0, typeO.cif.025), col=4, type='s', lwd=2, lty=2)
## points(c(0, post$times)*30.5, c(0, typeO.cif.975), col=4, type='s', lwd=2, lty=2)
##      legend(450, .4, c("Transplant(BART)", "Transplant(AJ)",
##                        "Death(AJ)", "Withdrawal(AJ)"),
##             col=c(4, 2, 1, 3), lwd=2)
