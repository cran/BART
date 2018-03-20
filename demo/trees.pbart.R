
library(BART)

f = function(x) #only the first 5 matter
    sin(pi*x[ , 1]*x[ , 2]) + 2*x[ , 3]*x[ , 4]^2 + x[ , 5]

N = 1000
sigma = 1.0  #y = f(x) + sigma*z where z~N(0, 1)
P = 10       #number of covariates

V = diag(P)
V[5, 6] = 0.8
V[6, 5] = 0.8
L <- chol(V)
set.seed(12)
x.train=matrix(rnorm(N*P), N, P) %*% L
dimnames(x.train)[[2]] <- paste0('x', 1:P)
round(cor(x.train), digits=2)

Ey.train = f(x.train)
y.train=((Ey.train+sigma*rnorm(N))>0)*1
table(y.train)

set.seed(21)
post = pbart(x.train, y.train, sparse=TRUE)
post$varprob.mean>1/P

##write(post$treedraws$trees, 'trees.pbart.txt')
tc <- textConnection(post$treedraws$tree)
trees <- read.table(file=tc, fill=TRUE,
                    row.names=NULL, header=FALSE,
                    col.names=c('node', 'var', 'cut', 'leaf'))
close(tc)
m <- 1 ## MCMC samples
M <- trees$node[1]
n <- 0 ## trees
H <- trees$var[1]
branch <- matrix(0, nrow=P, ncol=P)
dimnames(branch)[[1]] <- paste0('x', 1:P)
dimnames(branch)[[2]] <- paste0('x', 1:P)
L <- nrow(trees)
for(l in 2:L) {
    if(is.na(trees$leaf[l])) {
        n <- n+1
        if(n>H) {
            n <- 1
            m <- m+1
        }
        C <- trees$node[l] ## nodes in tree
        B <- (C-1)/2 ## branches in tree
        i <- 0
        j <- 0
        if(C>1) vars <- integer(C)
        branch. <- 0*branch
    }
    else if(B>1) {
        i <- i+1
        h <- trees$node[l]
        if(i<C) {
            t <- floor(log2(h))
            k <- h-2^t
            if(trees$node[l+1]==(2^(t+1)+2*k)) {
                vars[h] <- trees$var[l]+1
                j <- j+1
                if(j>B) stop('Too many branches')
            }
        }
        else {
            for(h. in (C-1):2) {
                h <- h.
                j <- vars[h]
                if(j!=0)
                    for(t in (floor(log2(h))-1):0) {
                        if((h%%2)==0) k <- (h-2^(t+1))/2
                        else k <- (h-2^(t+1)-1)/2
                        h <- 2^t+k
                        i <- vars[h]
                        if(i!=j) branch.[min(i, j), max(i, j)] <- 1
                        vars[h] <- 0
                    }
            }
            branch <- branch+branch.
        }
    }
}
C <- sum(c(branch))
for(i in 1:(P-1))
    for(j in (i+1):P)
        if(i!=j) branch[j, i] <- branch[i, j]/C
round(branch, digits=2)

