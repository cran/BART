
#include "common.h"
#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "heterbart.h"
#include "latent.h"
#include "rand_draws.h"

#ifndef NoRcpp

#define FHAT(a, b) fhat(a, b)

RcppExport SEXP cmbart(
   SEXP _ix,		//training data x, pxn
   SEXP _iy,		//training data y, nx1
   SEXP _iC,		//number of classes, y in {0,1,...C-1}
   SEXP _im,		//number of trees in BART sum
   SEXP _inc,		//number of cutpoints for decision rules
   SEXP _ind,           //number of kept draws (except for thinnning ..)
   SEXP _iburn,         //number of burn-in draws skipped
   SEXP _ipower,	//tree prior
   SEXP _ibase,		//tree prior
   SEXP _itau		//bottom node mu prior
)
{

   //--------------------------------------------------
   //process arguments
   Rcpp::NumericVector  xv(_ix);
   double *ix = &xv[0];
   Rcpp::IntegerVector yv(_iy);

   size_t m = Rcpp::as<int>(_im);
   size_t nc = Rcpp::as<int>(_inc);
   size_t nd = Rcpp::as<int>(_ind);
   size_t burn = Rcpp::as<int>(_iburn);
   double mybeta = Rcpp::as<double>(_ipower);
   double alpha = Rcpp::as<double>(_ibase);
   double tau = Rcpp::as<double>(_itau);

   size_t C = Rcpp::as<int>(_iC);

   //--------------------------------------------------
   //print out arguments (so you can check)
   size_t n = yv.size();
   size_t ntimesp = xv.size();
   size_t p=ntimesp/n;

   //return storage
   Rcpp::NumericMatrix fhat(n,C); //posterior mean of (f^j(x_i))
   Rcpp::NumericVector fdrawstrain = Rcpp::NumericVector(Rcpp::Dimension(nd,n,C));
   //to get (i,j,k) use k*d1*d2 + j*d1 + i

#else

#define FHAT(a, b) fhat[a][b]

void cmbart(
   double *ix,		//training data x, pxn
   int *yv,		//training data y, nx1
   int C,		//number of classes, y in {0,1,...C-1}
   size_t m,		//number of trees in BART sum
   size_t nc,		//number of cutpoints for decision rules
   size_t nd,           //number of kept draws (except for thinnning ..)
   size_t burn,         //number of burn-in draws skipped
   double mybeta,	//tree prior
   double alpha,        //tree prior
   double tau,		//bottom node mu prior
   size_t n,
   size_t p,
   double* _fhat,
   double* fdrawstrain
)
{

   std::vector<double*> fhat(n);

   for(size_t i=0; i<n; ++i) fhat[i]=&_fhat[i*C];

#endif

   printf("*****Into main of Multinomial BART\n");
   //-----------------------------------------------------------
   //random number generation
   //GetRNGstate();
   arn gen;
   newRNGstates();  /* for Bobbys OpenMP latents */

   cout << "n: " << n << endl;
//   cout << "n*p: " << ntimesp << endl;
   cout << "p: " << p << endl;
   cout << "Number classes: " << C << endl;

   printf("*****Number of Trees: %d\n",m);
   printf("*****Number of Cut Points: %d\n",nc);
   printf("*****burn and ndpost: %d, %d\n",burn,nd);
   printf("*****Prior:\nbeta,alpha,tau: %lf,%lf,%lf\n",mybeta,alpha,tau);

   //--------------------------------------------------


   //--------------------------------------------------
   //allocate latents and temporary storage
   double *z = new double[n*C]; //z plus or minus
   double *zpm = new double[n*C];
   double *lambda = new double[n*C];
   double *yf = new double[n];
   double *svec = new double[n];

   for(size_t i=0;i<C;i++) {
      for(size_t j=0;j<n;j++) {
         if(yv[j]==i) z[i*n+j]=1.0;
         else z[i*n+j]=-1.0;
      }
   }
   for(size_t i=0;i<(n*C);i++)  {
      lambda[i]=1.0; 
   }
   for(size_t i=0;i<n;i++) svec[i]=1.0;
   for(size_t i=0;i<(n*C);i++) zpm[i]=z[i];


   //--------------------------------------------------
   //set up barts
   std::vector<heterbart> bm(C);
   for(size_t i=0;i<C;i++) {
      bm[i].setm(m);
      bm[i].setprior(alpha,mybeta,tau);
      bm[i].setdata(p,n,ix,z+i*n,nc);
      bm[i].draw(svec,gen);
   }

   //--------------------------------------------------
   //mcmc loop
   size_t printevery=1;
   for(size_t i=0;i<(nd+burn);i++) {

      if(i%printevery==0) printf("done %d (out of %d)\n",i,nd+burn);

      for(size_t c=0;c<C;c++) {
         for(size_t j=0;j<n;j++) {
            yf[j]=zpm[c*n+j]*bm[c].f(j);
            double tempC=0.0;
            for(size_t k=0;k<C;k++) {
               if(k==c) continue;
               tempC+= exp(zpm[k*n+j]*bm[k].f(j));
            }
            yf[j] -= log(tempC);
         }
         draw_z(n,yf,lambda+n*c,z+n*c);
         for(size_t j=0;j<n;j++) z[c*n+j] *=zpm[c*n+j];
         draw_lambda(n,yf,1000,1.0,lambda+n*c);
         for(size_t j=0;j<n;j++) svec[j]=sqrt(lambda[n*c+j]);
         bm[c].draw(svec,gen);
      }
      if(i>=burn) {
         for(size_t c=0;c<C;c++) {
            for(size_t j=0;j<n;j++) {
               FHAT(j,c) += bm[c].f(j);
               fdrawstrain[c*nd*n+j*nd+i-burn]=bm[c].f(j);
            }
         }
      }
      
   }
   for(size_t c=0;c<C;c++) {
      for(size_t j=0;j<n;j++) FHAT(j,c) /= nd;
   }

   delete [] z;
   delete [] lambda;
   delete [] yf;
   delete [] svec;

   //PutRNGstate();
   deleteRNGstates();
   //--------------------------------------------------

#ifndef NoRcpp
   Rcpp::List ret;
   ret["fmeantrain"]=fhat;
   ret["fdrawstrain"]=fdrawstrain;
   return(ret);
#endif

}
