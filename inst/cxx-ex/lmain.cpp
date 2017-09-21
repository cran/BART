
#include "clbart.cpp"

int main(void) {
  size_t n=10, p=1, np=1, m=50, nc=100, nd=1000, burn=250,
    nkeeptrain=nd, nkeeptest=nd, nkeeptestmean=nd, nkeeptreedraws=nd,
    printevery=100;
  int treesaslists=0;
  unsigned int n1=111, n2=222;
  double k=2., power=2., mybeta=power, base=0.95, alpha=base, 
    binaryOffset=.0, tau=3./(k*sqrt(m));
  double xtrain[10]={0., 0., 0., 0., 0., 1., 1., 1., 1., 1.}, xtest[1]={0.};
  int y[10]={1, -1, -1, 1, -1, -1, 1, 1, -1, 1};
  double* trmean=new double[n];
  double* temean=new double[np];
//  double* sdraw=new double[nd+burn];
  double* _trdraw=new double[nkeeptrain*n];
  double* _tedraw=new double[nkeeptest*np];

   std::vector<double*> trdraw(nkeeptrain);
   std::vector<double*> tedraw(nkeeptest);

   for(size_t i=0; i<nkeeptrain; ++i) trdraw[i]=&_trdraw[i*n];
   for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*np];

  clbart(n, p, np, &xtrain[0], &y[0], &xtest[0], m, nc, nd, burn,
	 mybeta, alpha, binaryOffset, tau, nkeeptrain, nkeeptest,
	 nkeeptestmean, nkeeptreedraws, printevery, treesaslists,
	 n1, n2, trmean, temean, _trdraw, _tedraw);

#ifdef RNG_random
  cout << "RNG_random" << '\n';
#elif defined (RNG_Rmath)
  cout << "RNG_Rmath" << '\n';
#endif
  
  cout << ::plogis(temean[0], 0., 1., 1., 0.) << '\n';
  cout << 1.-::plogis(temean[0], 0., 1., 1., 0.) << '\n';
  //for(size_t i=0; i<np; ++i) cout << ::plogis(temean[i], 0., 1., 1., 0.) << '\n';

  for(size_t i=0; i<5; ++i) 
    for(size_t j=0; j<np; ++j) cout << i << '\t' << j << '\t' << ::plogis(tedraw[i][j], 0., 1., 1., 0.) << '\n';

  delete[] trmean;
  delete[] temean;
//  delete[] sdraw;
  delete[] _trdraw;
  delete[] _tedraw;

  return 0;
}
