
#include "clbart.cpp"

int main(void) {
  size_t n=10, p=1, np=2, m=50, nd=1000, burn=250,
    nkeeptrain=nd, nkeeptest=nd, nkeeptreedraws=nd,
    printevery=100;
  unsigned int n1=111, n2=222;
  double k=2., power=2., mybeta=power, base=0.95, alpha=base, 
    binaryOffset=.0, tau=3./(k*sqrt(m));
  double xtrain[10]={0., 0., 0., 0., 0., 1., 1., 1., 1., 1.}, xtest[2]={0., 1.};
  int y[10]={1, -1, -1, 1, -1, -1, 1, 1, -1, 1};
  double* _trdraw=new double[nkeeptrain*n];
  double* _tedraw=new double[nkeeptest*np];
  int nc[1]={100};

   std::vector<double*> trdraw(nkeeptrain);
   std::vector<double*> tedraw(nkeeptest);

   for(size_t i=0; i<nkeeptrain; ++i) trdraw[i]=&_trdraw[i*n];
   for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*np];

  clbart(n, p, np, &xtrain[0], &y[0], &xtest[0], m, &nc[0], nd, burn,
	 mybeta, alpha, binaryOffset, tau, //0, 0., 0., NULL, 0., 0., 0., 0, 
	 0, 0., 0., 0., 0, 
	 nkeeptrain, nkeeptest, nkeeptreedraws, printevery, 
	 n1, n2, _trdraw, _tedraw);

#ifdef RNG_random
  cout << "RNG_random" << '\n';
#elif defined (RNG_Rmath)
  cout << "RNG_Rmath" << '\n';
#endif
  
  double prob_mean0=0., prob_mean1=0.;

  for(size_t i=0; i<nd; ++i) {
    prob_mean0 += ::plogis(tedraw[i][0], 0., 1., 1., 0.)/nd;
    prob_mean1 += ::plogis(tedraw[i][1], 0., 1., 1., 0.)/nd;
  }

  cout << "P(y=1, x=0)=" << prob_mean0 << '\n';
  cout << "P(y=1, x=1)=" << prob_mean1 << '\n';

  delete[] _trdraw;
  delete[] _tedraw;

  return 0;
}
