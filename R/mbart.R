mbart = function(
   x.train,			#training data x, nxp 
   y.train,			#training data y, nx1
   k=2.0,
   power=2.0, 
   base=.95,
   ntree=200L,
   numcut=100L,
   ndpost=1000L, 
   nskip=100L
)
{
#--------------------------------------------------
tau=(3.0)/(2*k*sqrt(ntree))
C = length(unique(y.train))
yy = as.integer(as.factor(y.train))-1 	#make it 0,1,,,,C-1
#--------------------------------------------------
res = .Call("cmbart",
   t(x.train),
   yy,
   C,
   ntree,
   numcut,
   ndpost,
   nskip,
   power,
   base,
   tau
)
attr(res, 'class')='mbart'
return(res)
}
