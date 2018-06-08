
#include "cwbart.cpp"

int main(void) {
  size_t n=10, p=1, np=2, m=50, nd=1000, burn=250,
    nkeeptrain=nd, nkeeptest=nd, nkeeptestmean=nd, nkeeptreedraws=nd,
    printevery=100;
  unsigned int n1=111, n2=222;
  double k=2., power=2., mybeta=power, base=0.95, alpha=base, 
    tau=14./(2.*k*sqrt(m)), nu=3., ymu=0., ysd=0., sigquant=0.9, lambda;
  double xtrain[10]={0., 0., 0., 0., 0., 1., 1., 1., 1., 1.}, xtest[2]={0., 1.};
  double y[10]={0., 1., -1., 2., -2., 10., 11., 9., 12., 8.};
  double w[10]={1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
  int nc[1]={100};

  for(size_t i=0; i<n; ++i) ymu += y[i]/n;
  for(size_t i=0; i<n; ++i) ysd += pow(y[i]-ymu, 2.)/(n-1.);
  lambda = ysd* ::qchisq(1.0-sigquant, nu, 1, 0)/nu;
  ysd=sqrt(ysd);

  double* trmean=new double[n];
  double* temean=new double[np];
  double* sdraw=new double[nd+burn];
  double* _trdraw=new double[nkeeptrain*n];
  double* _tedraw=new double[nkeeptest*np];

   std::vector<double*> trdraw(nkeeptrain);
   std::vector<double*> tedraw(nkeeptest);

   for(size_t i=0; i<nkeeptrain; ++i) trdraw[i]=&_trdraw[i*n];
   for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*np];

  cwbart(n, p, np, &xtrain[0], &y[0], &xtest[0], m, &nc[0], nd, burn,
	 mybeta, alpha, tau, nu, lambda, ysd, w, 0, 0., 0., NULL, 0., 0., 0., 0,
	 nkeeptrain, nkeeptest, nkeeptestmean, nkeeptreedraws, printevery, 
	 n1, n2, trmean, temean, sdraw, _trdraw, _tedraw);

#ifdef RNG_random
  cout << "RNG_random" << '\n';
#elif defined (RNG_Rmath)
  cout << "RNG_Rmath" << '\n';
#endif
  
  cout << temean[0] << '\n';
  cout << temean[1] << '\n';

//  for(size_t i=0; i<5; ++i) 
//    for(size_t j=0; j<np; ++j) cout << i << '\t' << j << '\t' << tedraw[i][j] << '\n';

  delete[] trmean;
  delete[] temean;
  delete[] sdraw;
  delete[] _trdraw;
  delete[] _tedraw;

  return 0;
}
