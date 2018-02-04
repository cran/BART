library(BART)

times <- matrix(c(3,  8,  9,
                  4, 12, 12),
                nrow=2, ncol=3, byrow=TRUE)

tstop <- matrix(c(7, 8, 0,
                  7, 0, 0),
                nrow=2, ncol=3, byrow=TRUE)

delta <- matrix(c(1, 1, 0,
                  1, 0, 0),
                nrow=2, ncol=3, byrow=TRUE)

data. <- recur.pre.bart(times=times, delta=delta,
                        tstop=tstop)
print(data.)

data. <- recur.pre.bart(times=times, delta=delta,
                        tstop=tstop, last.value=FALSE)
print(data.$tx.test)
