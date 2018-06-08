/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017-2018 Robert McCulloch and Rodney Sparapani
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

#include "common.h"
#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
/*
#include "heterbart.h"
#include "latent.h"
#include "rand_draws.h"
*/
#include "rtnorm.h"

#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)
#define VARCNT(a, b) varcnt(a, b)
#define VARPRB(a, b) varprb(a, b)

RcppExport SEXP cmbart(
   SEXP _in,            //number of observations in training data
   SEXP _ip,            //dimension of x
   SEXP _inp,           //number of observations in test data
   SEXP _ix,            //x, train,  pxn (transposed so rows are contiguous in memory)
   SEXP _iy,            //y, train,  nx1
   SEXP _iC,		//number of classes, y in {0,1,...C-1}
   SEXP _ixp,           //x, test, pxnp (transposed so rows are contiguous in memory)
   SEXP _im,            //number of trees
   SEXP _inc,           //number of cut points
   SEXP _ind,           //number of kept draws (except for thinnning ..)
   SEXP _iburn,         //number of burn-in draws skipped
   SEXP _ipower,
   SEXP _ibase,
   SEXP _binaryOffset,
   SEXP _itau,
   SEXP _idart,         //dart prior: true(1)=yes, false(0)=no
   SEXP _itheta,
   SEXP _iomega,
   SEXP _igrp,
   SEXP _ia,            //param a for sparsity prior
   SEXP _ib,            //param b for sparsity prior
   SEXP _irho,          //param rho for sparsity prior (default to p)
   SEXP _iaug,          //categorical strategy: true(1)=data augment false(0)=degenerate trees
   SEXP _inkeeptrain,
   SEXP _inkeeptest,
//   SEXP _inkeeptestme,
   SEXP _inkeeptreedraws,
   SEXP _inprintevery,
//   SEXP _keepevery,
   //SEXP _treesaslists,
   SEXP _Xinfo
)
{

   //--------------------------------------------------
   //process arguments
   size_t n = Rcpp::as<int>(_in);
   size_t p = Rcpp::as<int>(_ip);
   size_t np = Rcpp::as<int>(_inp);
   Rcpp::NumericVector  xv(_ix);
   double *ix = &xv[0];
   Rcpp::IntegerVector  yv(_iy); // O, ..., C-1
   int *iy = &yv[0];
   Rcpp::NumericVector  xpv(_ixp);
   double *ixp = &xpv[0];
   size_t m = Rcpp::as<int>(_im);
   //size_t nc = Rcpp::as<int>(_inc);
   Rcpp::IntegerVector _nc(_inc);
   int *numcut = &_nc[0];
   size_t nd = Rcpp::as<int>(_ind);
   size_t burn = Rcpp::as<int>(_iburn);
   double mybeta = Rcpp::as<double>(_ipower);
   double alpha = Rcpp::as<double>(_ibase);

   // mbart does not currently employ the binaryOffset trick
   //double binaryOffset = Rcpp::as<double>(_binaryOffset);
   Rcpp::NumericVector bO(_binaryOffset);
   double *binaryOffset = &bO[0];

   double tau = Rcpp::as<double>(_itau);
   bool dart;
   if(Rcpp::as<int>(_idart)==1) dart=true;
   else dart=false;
   double a = Rcpp::as<double>(_ia);
   double b = Rcpp::as<double>(_ib);
   double rho = Rcpp::as<double>(_irho);
   bool aug;
   if(Rcpp::as<int>(_iaug)==1) aug=true;
   else aug=false;
   double theta = Rcpp::as<double>(_itheta);
   double omega = Rcpp::as<double>(_iomega);
   size_t nkeeptrain = Rcpp::as<int>(_inkeeptrain);
   size_t nkeeptest = Rcpp::as<int>(_inkeeptest);
//   size_t nkeeptestme = Rcpp::as<int>(_inkeeptestme);
   size_t nkeeptreedraws = Rcpp::as<int>(_inkeeptreedraws);
   size_t printevery = Rcpp::as<int>(_inprintevery);
//   int keepevery = Rcpp::as<int>(_keepevery);
   //int treesaslists = Rcpp::as<int>(_treesaslists);
   Rcpp::NumericMatrix Xinfo(_Xinfo);
//   Rcpp::IntegerMatrix varcount(nkeeptreedraws, p);
   size_t C = Rcpp::as<int>(_iC);
   Rcpp::NumericMatrix varprb(nkeeptreedraws,C*p);
   Rcpp::IntegerMatrix varcnt(nkeeptreedraws,C*p);

   //return storage
   Rcpp::NumericMatrix trdraw(nkeeptrain,n*C), tedraw(nkeeptest,np*C); 
   //Rcpp::NumericMatrix trmean(n,C); //posterior mean of (f^j(x_i))
   //Rcpp::NumericVector trdraw = Rcpp::NumericVector(Rcpp::Dimension(nd,n,C));
   //to get (i,j,k) use k*d1*d2 + j*d1 + i

   xinfo _xi;
   if(Xinfo.size()>0) {
     _xi.resize(p);
     for(size_t i=0;i<p;i++) {
       _xi[i].resize(numcut[i]);
       for(size_t j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
     }
   }

   arn gen;

#else

#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]
#define VARCNT(a, b) varcnt[a][b]
#define VARPRB(a, b) varprb[a][b]

void cmbart(
   size_t n,            //number of observations in training data
   size_t p,		//dimension of x
   size_t np,		//number of observations in test data
   double* ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
   int* iy,		//y, train,  nx1
   size_t C,		//number of classes, y in {0,1,...C-1}
   double* ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
   size_t m,		//number of trees
   int *numcut,		//number of cut points
   size_t nd,		//number of kept draws (except for thinnning ..)
   size_t burn,		//number of burn-in draws skipped
   double mybeta,
   double alpha,
   double binaryOffset,
   double tau,
   bool dart,           //dart prior: true(1)=yes, false(0)=no        
   double theta,
   double omega,
   int *grp,
   double a,		//param a for sparsity prior                 
   double b,		//param b for sparsity prior                 
   double rho,		//param rho for sparsity prior (default to p) 
   bool aug,		//categorical strategy: true(1)=data augment false(0)=degenerate trees
   size_t nkeeptrain,
   size_t nkeeptest,
//   size_t nkeeptestme,
   size_t nkeeptreedraws,
   size_t printevery,
//   int keepevery,
   //int treesaslists,
   unsigned int n1, // additional parameters needed to call from C++
   unsigned int n2,
   double* _trdraw,
   double* _tedraw
)
{

  double** trdraw=new double*[nkeeptrain];
  for(size_t i=0; i<nkeeptrain; ++i) trdraw[i]=&_trdraw[i*C*n];

  double** tedraw=new double*[nkeeptest];
  for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*C*np];

//  double** trmean=new double*[n];
//  for(size_t i=0; i<n; ++i) trmean[i]=&_trmean[i*C];

   //matrix to return dart posteriors (counts and probs)
   std::vector< std::vector<size_t> > varcnt;
   std::vector< std::vector<double> > varprb;

   arn gen(n1, n2);
#endif

   double* fhattest=0; //out-of-sample fit
   if(np) { fhattest = new double[np]; }
   else nkeeptest=0;

   size_t skiptr,skipte,/*skipteme,*/skiptreedraws;
   if(nkeeptrain) {skiptr=nd/nkeeptrain;}
   else skiptr = nd+1;
   if(nkeeptest) {skipte=nd/nkeeptest;}
   else skipte=nd+1;
//   if(nkeeptestme) {skipteme=nd/nkeeptestme;}
//   else skipteme=nd+1;
   if(nkeeptreedraws) {skiptreedraws = nd/nkeeptreedraws;}
   else skiptreedraws=nd+1;

//   size_t trcnt=0; //count kept train draws
//   size_t tecnt=0; //count kept test draws
//   size_t temecnt=0; //count test draws into posterior mean
//   size_t treedrawscnt=0; //count kept bart draws
   bool keeptest,keeptrain,/*keeptestme*/keeptreedraw;

   printf("*****Into main of Multinomial BART\n");
   //-----------------------------------------------------------
   //random number generation
   //newRNGstates();  /* for Bobbys OpenMP latents */
/*
   cout << "n: " << n << endl;
//   cout << "n*p: " << ntimesp << endl;
   cout << "p: " << p << endl;
   cout << "Number classes: " << C << endl;
*/
   printf("*****Number of Trees: %zu\n",m);
   printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p-1]);
   printf("*****burn and ndpost: %zu, %zu\n",burn,nd);
   printf("*****Prior:mybeta,alpha,tau: %lf,%lf,%lf\n",
                   mybeta,alpha,tau);
   printf("*****binaryOffset: %lf\n",binaryOffset);
   cout << "*****Dirichlet:sparse,theta,omega,a,b,rho,augment: " 
	<< dart << ',' << theta << ',' << omega << ',' << a << ',' 
	<< b << ',' << rho << ',' << aug << endl;

   //--------------------------------------------------
   //allocate latents and temporary storage
   double *z = new double[n]; //z plus or minus
//   double *z = new double[n*C]; //z plus or minus
/*
   double *zpm = new double[n*C];
   double *lambda = new double[n*C];
   double *yf = new double[n];
   double *svec = new double[n];
   for(size_t i=0;i<C;i++) {
      for(size_t j=0;j<n;j++) {
         if(iy[j]==i) z[i*n+j]=1.0;
         else z[i*n+j]=-1.0;
      }
   }
   for(size_t i=0;i<(n*C);i++)  {
      lambda[i]=1.0; 
   }
   for(size_t i=0;i<n;i++) svec[i]=1.0;
   for(size_t i=0;i<(n*C);i++) zpm[i]=z[i];
*/


   //--------------------------------------------------
   //set up barts
   std::vector<bart> bm(C);
   std::vector<std::stringstream> trees(C);  //string stream to write trees to
          
   for(size_t i=0;i<C;i++) {
     trees[i].precision(10);
     trees[i] << nkeeptreedraws << " " << m << " " << p << endl;
      bm[i].setm(m);
      bm[i].setprior(alpha,mybeta,tau);
      bm[i].setdata(p,n,ix,z,numcut);
      //bm[i].setdata(p,n,ix,z+i*n,numcut);
      bm[i].setdart(a,b,rho,aug,dart);
#ifndef NoRcpp
      if(Xinfo.size()>0) bm[i].setxinfo(_xi);
#endif
      for(size_t k=0; k<n; k++) {
	if(iy[k]!=i) z[k]= -rtnorm(0., binaryOffset[i], 1., gen);
	else z[k]=rtnorm(0., -binaryOffset[i], 1., gen);
      }

      bm[i].draw(1., gen);
   }

   // dart iterations
   std::vector<double> ivarprb (p,0.);
   std::vector<size_t> ivarcnt (p,0);

   xinfo& xi = bm[0].getxinfo(); // they should all be alike

   //--------------------------------------------------
   //mcmc loop
   for(size_t i=0;i<(nd+burn);i++) {

      if(i%printevery==0) printf("done %zu (out of %zu)\n",i,nd+burn);

      for(size_t c=0;c<C;c++) {
	if(i==(burn/2)&&dart) bm[c].startdart();

         for(size_t k=0; k<n; k++) {
	   if(iy[k]!=c) z[k]= -rtnorm(-bm[c].f(k), binaryOffset[c], 1., gen);
	   else z[k]=rtnorm(bm[c].f(k), -binaryOffset[c], 1., gen);
	 }

         bm[c].draw(1., gen);
      }
/*
      for(size_t c=0;c<C;c++) {
	if(i==(burn/2)&&dart) bm[c].startdart();
         for(size_t j=0;j<n;j++) {
            yf[j]=zpm[c*n+j]*bm[c].f(j);
            double tempC=0.0;
            for(size_t k=0;k<C;k++) {
               if(k==c) continue;
               tempC+= exp(zpm[k*n+j]*bm[k].f(j));
            }
            //yf[j] -= log(tempC);
         }
         draw_z(n,yf,lambda+n*c,z+n*c);
         //for(size_t j=0;j<n;j++) z[c*n+j] *=zpm[c*n+j];
         draw_lambda(n,yf,1000,1.0,lambda+n*c);
         //for(size_t j=0;j<n;j++) svec[j]=sqrt(lambda[n*c+j]);
         for(size_t j=0;j<n;j++) {
	   z[c*n+j] *=zpm[c*n+j];
	   svec[j]=sqrt(lambda[n*c+j]);
	 }
         bm[c].draw(svec,gen);
      }
*/

      keeptest = nkeeptest && (((i-burn+1) % skipte) ==0);
      keeptrain = nkeeptrain && (((i-burn+1) % skiptr) ==0);
      keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);

      if(i>=burn && (keeptest || keeptrain || keeptreedraw)) {
	for(size_t c=0;c<C;c++) {
	  if(keeptreedraw) {
	    for(size_t j=0;j<m;j++) trees[c] << bm[c].gettree(j);
	    ivarcnt=bm[c].getnv();
	    ivarprb=bm[c].getpv();
	    size_t k=(i-burn)/skiptreedraws;
	    for(size_t j=0;j<p;j++){
	      VARCNT(k,j*C+c)=ivarcnt[j];
	      VARPRB(k,j*C+c)=ivarprb[j];
	    }
	  }
	  if(keeptrain) {
	    size_t k=(i-burn)/skiptr;
	    for(size_t j=0;j<n;j++) 
	      TRDRAW(k, j*C+c)=binaryOffset[c]+bm[c].f(j);
	  }
	  //trdraw[c*nd*n+j*nd+i-burn]=bm[c].f(j);
	  //TRMEAN(j,c) += bm[c].f(j)/nd;
	  if(keeptest) {
	    size_t k=(i-burn)/skipte;
	    bm[c].predict(p,np,ixp,fhattest);
	    for(size_t j=0;j<np;j++) 
	      TEDRAW(k, j*C+c)=binaryOffset[c]+fhattest[j];
	  }
	}
      }
   }

   if(fhattest) delete [] fhattest;
   delete [] z;
   //delete [] lambda;
   //delete [] yf;
   //delete [] svec;

   //deleteRNGstates();
   //--------------------------------------------------

#ifndef NoRcpp
   Rcpp::List ret;
   ret["yhat.train"]=trdraw;
   if(np) ret["yhat.test"]=tedraw;
   ret["varcount"]=varcnt;
   ret["varprob"]=varprb;

   Rcpp::List xiret(xi.size());
   for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
   }

   Rcpp::List treesL;
   treesL["cutpoints"] = xiret;
   for(size_t c=0; c<C; ++c) {
     char name[6]="trees";
     name[4]=49+c;
     treesL[name]=Rcpp::CharacterVector(trees[c].str());
   }
   //if(treesaslists) treesL["lists"]=list_of_lists;
   ret["treedraws"] = treesL;

   return(ret);
#else
   delete [] trdraw;
   delete [] tedraw;
#endif

}
