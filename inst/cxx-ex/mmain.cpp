
#include "cmbart.cpp"

int main(void) {
  size_t C=2, n=10, p=1, np=2, m=50, nd=1000, burn=250,
    nkeeptrain=nd, nkeeptest=nd, nkeeptestmean=nd, nkeeptreedraws=nd,
    printevery=100;
  int treesaslists=0;
  unsigned int n1=111, n2=222;
  double k=2., power=2., mybeta=power, base=0.95, alpha=base, 
    binaryOffset=.0, tau=3./(k*sqrt(m));
  double xtrain[10]={0., 0., 0., 0., 0., 1., 1., 1., 1., 1.}, xtest[2]={0., 1.};
  int y[10]={1, 0, 0, 1, 0, 0, 1, 1, 0, 1};
  int nc[1]={100};

  double* _trdraw=new double[nkeeptrain*C*n];
  double* _tedraw=new double[nkeeptest*C*np];

   std::vector<double*> trdraw(nkeeptrain);
   std::vector<double*> tedraw(nkeeptest);

   for(size_t i=0; i<nkeeptrain; ++i) trdraw[i]=&_trdraw[i*C*n];
   for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*C*np];

  cmbart(n, p, np, &xtrain[0], &y[0], C, &xtest[0], m, &nc[0], nd, burn,
	 mybeta, alpha, binaryOffset, tau, false, 1., 1., p, true, 
	 nkeeptrain, nkeeptest, nkeeptreedraws, printevery, 
	 n1, n2, _trdraw, _tedraw);

  delete[] _trdraw;
  delete[] _tedraw;

  return 0;
}
