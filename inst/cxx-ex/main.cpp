
#include "common.h"
#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "heterbart.h"

int main(void) {
  /* 
     N: size of training data
     M: size of testing data
     P: covariates
     H: trees
     D: MCMC draws
     B: burn-in discarded
     C: cutpoints (an array since some variables may be discrete)
  */

  size_t N=10, M=2, P=1, H=50, B=0, D=1, T=B+D; // B=250, D=1000,
  int C[1]={100}, augment=0, sparse=0, theta=0, omega=1;
  unsigned int n1=111, n2=222; // seeds
  double k=2., power=2., base=0.95, a=0.5, b=1., tau=14./(2.*k*sqrt(H)), 
    nu=3., ymu=0., ysd=0., sigquant=0.9, lambda, sigma=1.;
  double xtrain[10]={0., 0., 0., 0., 0., 1., 1., 1., 1., 1.}, 
    xtest[2]={0., 1.};
  double y[10]={0., 1., -1., 2., -2., 10., 11., 9., 12., 8.}; // outcome
  double w[10]={1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}; // known weights for SDi

  // tune lambda parameter 
  for(size_t i=0; i<N; ++i) ymu += y[i]/N;
  for(size_t i=0; i<N; ++i) ysd += pow(y[i]-ymu, 2.)/(N-1.);
  lambda = ysd* ::qchisq(1.-sigquant, nu, 1, 0)/nu;
  ysd=sqrt(ysd);

  double* sd=new double[B+D];
  double* _ftrain=new double[N*D];
  double* _ftest =new double[M*D];

  std::vector<double*> ftrain(D);
  std::vector<double*> ftest(D);

  for(size_t h=0; h<D; ++h) ftrain[h]=&_ftrain[h*N];
  for(size_t h=0; h<D; ++h) ftest[h] =&_ftest[h*M];

  arn gen(n1, n2); // abstract random number generator

  heterbart A(H);
  A.setprior(base, power, tau);
  A.setdata(P, N, xtrain, y, C);
  A.setdart(a, b, P, augment, sparse, theta, omega);

  double *svec = new double[N]; //sigma vector
  for(size_t i=0; i<N; i++) svec[i]=w[i]*sigma;

  double* fhattest=0; //prediction
  if(M) fhattest = new double[M]; 

  for(size_t h=0; h<T; ++h) {
    //if(h==(burn/2)&&dart) A.startdart();
    A.draw(svec, gen); //draw f

    //draw sigma
    double restemp, rss;
    rss=0.;
    for(size_t i=0; i<N; i++) {
      restemp=(y[i]-A.f(i))/w[i]; 
      rss += restemp*restemp;
    }
    sigma = sqrt((nu*lambda + rss)/gen.chi_square(N+nu));
    sd[h]=sigma;

    for(size_t i=0; i<N; i++) svec[i]=w[i]*sigma;

    size_t j;
    j=h-B;

    if(j>=0) {
      for(size_t i=0; i<N; i++) ftrain[j][i]=A.f(i);
      A.predict(P, M, xtest, fhattest); //predict
      for(size_t i=0; i<M; i++) ftest[j][i]=fhattest[i];
    }
  }

  cout << "Draw 1 MCMC sample from the BART prior\n";

  cout << "f(x) train \n";
  for(size_t i=0; i<N; ++i) 
    cout << "i=" << i << ", x=" << xtrain[i] << ", y=" << y[i] << ",\tf(x)=" << ftrain[0][i] << '\n';

  cout << "f(x) test \n";
  for(size_t i=0; i<M; ++i) 
    cout << "i=" << i << ", x=" << xtest[i] << ", f(x)=" << ftest[0][i] << '\n';

  delete[] sd;
  delete[] _ftrain;
  delete[] _ftest;
  delete [] svec;
  if(fhattest) delete[] fhattest;

  return 0;
}
