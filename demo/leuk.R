
library(BART)

data(leukemia)

leukemia$TD=ceiling(leukemia$TD/30)
leukemia$TB=ceiling(leukemia$TB/30)
leukemia$TA=ceiling(leukemia$TA/30)
leukemia$TC=ceiling(leukemia$TC/30)
leukemia$TP=ceiling(leukemia$TP/30)
leukemia$X7=ceiling(leukemia$X7/30)

N=137
events=unique(sort(c(leukemia$TD, leukemia$TB)))
K=length(events)
T=events[K]
    ## the following covariates appear to be important
    ## G, TB, R, TA, A, TC,  C, TP,  P, X2, X8
pick=c(1,  3, 5,  7, 8,  9, 10, 11, 12, 14, 20)
x.train3=as.matrix(leukemia[ , -c(2, 4, 6)])
P=ncol(x.train3)
L=32
x.test3=matrix(nrow=L*N, ncol=P)
dimnames(x.test3)=dimnames(x.train3)
k=1
for(R in 0:1)
    for(A in 0:1)
        for(C in 0:1)
            for(P in 0:1)
                for(X2 in c(23, 32)) {
                    h=(k-1)*N+1:N
                    x.test3[h, ]=x.train3
                    x.test3[h, 'TB']=R*8+(1-R)*T
                    x.test3[h, 'R']=R
                    x.test3[h, 'TA']=A*1+(1-A)*T
                    x.test3[h, 'A']=A
                    x.test3[h, 'TC']=C*5+(1-C)*T
                    x.test3[h, 'C']=C
                    x.test3[h, 'TP']=P*1+(1-P)*T
                    x.test3[h, 'P']=P
                    x.test3[h, 'X2']=X2
                    k=k+1
                }

post3=mc.surv.bart(x.train=x.train3, times=leukemia$TD, delta=leukemia$D, 
                   events=events, ztimes=c(2, 4, 6, 8), zdelta=c(3, 5, 7, 9),
                   sparse=TRUE, mc.cores=8, seed=99)

state3=surv.pre.bart(leukemia$TD, leukemia$D, x.train3, x.test3,
                     events=events,
                     ztimes=c(2, 4, 6, 8), zdelta=c(3, 5, 7, 9))

## post3=mc.surv.bart(state3$tx.train, state3$y.train,
##                    x.test=state3$tx.train,
##                    sparse=TRUE,
##                    mc.cores=8, seed=99)

x.train2=x.train3[ , -(2:3)]
x.test2=x.test3[ , -(2:3)]
post2=mc.surv.bart(x.train=x.train2, times=leukemia$TB, delta=leukemia$R, 
                   events=events, ztimes=c(2, 4, 6), zdelta=c(3, 5, 7),
                   sparse=TRUE, mc.cores=8, seed=99)

state2=surv.pre.bart(leukemia$TB, leukemia$R, x.train2, x.test2,
                     events=events, ztimes=c(2, 4, 6), zdelta=c(3, 5, 7))

## post2=mc.surv.bart(state2$tx.train, state2$y.train,
##                    x.test=state2$tx.train,
##                    sparse=TRUE,
##                    mc.cores=8, seed=99)

##pdf(file='leuk.pdf')
par(mfrow=c(2, 2))
for(l in 1:L) {
    h=(l-1)*N*K+1:(N*K)

    for(G in 1:5) {
        if(G==1) {
            state3$tx.test[h, 'G']=1
            state3$tx.test[h, 'X8']=0
        } else if(G %in% 2:3) {
            state3$tx.test[h, 'G']=2
            state3$tx.test[h, 'X8']=G-2
        } else if(G %in% 4:5) {
            state3$tx.test[h, 'G']=3
            state3$tx.test[h, 'X8']=G-4
        } 
        state2$tx.test[h, 'G']=state3$tx.test[h, 'G']
        state2$tx.test[h, 'X8']=state3$tx.test[h, 'X8']
        pred3=predict(post3, state3$tx.test[h, ], mc.cores=8)
        pred2=predict(post2, state2$tx.test[h, ], mc.cores=8)
        i=(l-1)*N+1
        R=x.test3[i, 'R']
        string=paste0(' R=', R, ' A=', x.test3[i, 'A'], 
                      ' C=', x.test3[i, 'C'], ' P=', x.test3[i, 'P'], 
                      ' X2=', x.test3[i, 'X2'])
        
        state0.mean=double(K)
        state1.mean=double(K)
        for(j in 1:K) {
            k=seq(j, N*K, by=K) 
            state0.mean[j]=mean(apply(pred3$surv.test[ , k], 1, mean))
            state1.mean[j]=mean(apply(pred2$surv.test[ , k]*
                                      pred3$surv.test[ , k], 1, mean))
            ## state2.mean[j]=mean(apply((1-pred2$surv.test[ , k])*
            ##                           pred3$surv.test[ , k], 1, mean))
        }

        if(R==1) state1.mean[8:K]=0
        
        if(G==1) 
        plot(c(0, pred3$times), c(1, state0.mean), type='s', lwd=2, lty=G,
             ylim=0:1, xlab='t (months)', ylab='P(t, state)', main=string)
        else
        lines(c(0, pred3$times), c(1, state0.mean), type='s', lwd=2, lty=G)
        lines(c(0, pred3$times), c(1, state1.mean), lty=G,
              type='s', lwd=2, col=2)
        lines(c(0, pred3$times), c(0, state0.mean-state1.mean), lty=G,
              type='s', lwd=2, col=4)
        if((l%%4)==0) {
            legend('topright', col=c(1, 2, 4), lty=1, lwd=2,
                   legend=c('alive', 'remission', 'relapsed'))
        }
        else if((l%%2)==0) {
            legend('topright', lty=1:5, lwd=2,
                   legend=c('G=1 X8=0', 'G=2 X8=0', 'G=2 X8=1',
                            'G=3 X8=0', 'G=3 X8=1')) 
        }
    }
}
dev.off()
##par(mfrow=c(1, 1))
