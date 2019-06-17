
#include "cdpmwbart.cpp"

int main(void) {
  size_t n=10, p=1, np=1, m=50, nc=100, nd=1000, burn=250, nkeeptrain=nd;
  unsigned int n1=111, n2=222;
  double k=2., power=2., mybeta=power, base=0.95, alpha=base, initm=0., inits=1.,
   range=6., tau=range/(2.*k*sqrt(m));
  double xtrain[10]={0., 0., 0., 0., 0., 1., 1., 1., 1., 1.}, xtest[1]={0.};
  double y[10]={-3., -2., -1., 0., 0., 0., 0., 1., 2., 3.};
  double w[10]={1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
  double nu=3., qchi=::qchisq(0.1, nu, 1, 0), lambda=pow(inits, 2.)*qchi/nu;

  std::vector<double> gm={-0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1.};
  std::vector<double> pm={0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
  std::vector<double> gs={0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2};
  std::vector<double> ps={0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

  double* sigmadr=new double[nd+burn];
  double* _trdraw=new double[nkeeptrain*n];
  double* trmean=new double[n];
  double* _mdraw=new double[nkeeptrain*n];
  double* _sdraw=new double[nkeeptrain*n];
  int* istar=new int[nkeeptrain];

  cdpmwbart(n, p, np, &xtrain[0], &y[0], &xtest[0], nd, burn, m, nc, mybeta, alpha, tau, initm, inits,
	   gm, 0.1, pm, gs, 0.1, ps, lambda, nu, w, 1, n1, n2,
	   sigmadr, _trdraw, trmean, _mdraw, _sdraw, istar);

#ifdef RNG_random
  cout << "RNG_random" << '\n';
#elif defined (RNG_Rmath)
  cout << "RNG_Rmath" << '\n';
#endif

  delete[] sigmadr;
  delete[] _trdraw;
  delete[] trmean;
  delete[] _mdraw;
  delete[] _sdraw;
  delete[] istar;

  return 0;
}
