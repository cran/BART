/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
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

#include <ctime>
#include "dps.h"
#include "dpm.h"

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "heterbart.h"

#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define MDRAW(a, b) mdraw(a, b)
#define SDRAW(a, b) sdraw(a, b)

RcppExport SEXP cdpmwbart(
   SEXP _in,            //number of observations in training data
   SEXP _ip,		//dimension of x
   SEXP _inp,		//number of observations in test data
   SEXP _x,
   SEXP _y,
   SEXP _xp,
   SEXP _nd,
   SEXP _burn,
   SEXP _m,
   SEXP _nc,
   SEXP _power,
   SEXP _base,
   SEXP _tau,
   SEXP _mstart,
   SEXP _sstart,
   SEXP _gm,
   SEXP _am,
   SEXP _pm,
   SEXP _gs,
   SEXP _as,
   SEXP _ps,
   SEXP _lambda,
   SEXP _nu,
   SEXP _drawmu,
   SEXP _drawsd,
   SEXP _iw,
   SEXP _iq
)
{

   arn gen;

   size_t n = Rcpp::as<int>(_in);
   size_t p = Rcpp::as<int>(_ip);
   size_t np = Rcpp::as<int>(_inp);

   Rcpp::NumericMatrix xm(_x);
   double *x = &xm[0];

   Rcpp::NumericVector yv(_y);
   double *y = &yv[0];

   Rcpp::NumericMatrix xpm(_xp);
   //double *xp = &xpm[0];

   Rcpp::NumericVector wv(_iw);
   double *iw = &wv[0];
   double iq = Rcpp::as<double>(_iq);

   size_t nd = Rcpp::as<int>(_nd);
   size_t burn = Rcpp::as<int>(_burn);
   size_t m = Rcpp::as<int>(_m);
   size_t nc = Rcpp::as<int>(_nc);

   double mybeta = Rcpp::as<double>(_power);
   double alpha = Rcpp::as<double>(_base);
   double tau = Rcpp::as<double>(_tau);

   double mstart = Rcpp::as<double>(_mstart);
   double sstart = Rcpp::as<double>(_sstart);

   //dp prior for mu_i
   Rcpp::NumericVector gm(_gm);
   Rcpp::NumericVector pm(_pm);
   size_t mgsize = gm.size();
   double am = Rcpp::as<double>(_am);

   //dp prior for sigma_i
   Rcpp::NumericVector gs(_gs);
   Rcpp::NumericVector ps(_ps);
   size_t sgsize = gs.size();
   double as = Rcpp::as<double>(_as);

   double lambda = Rcpp::as<double>(_lambda);
   double nu = Rcpp::as<double>(_nu);

   //flags for what should be drawn
   size_t drawmu = Rcpp::as<int>(_drawmu);
   size_t drawsd = Rcpp::as<int>(_drawsd);

   Rcpp::NumericVector trmean(n); //train
   Rcpp::NumericMatrix trdraw(nd,n);
   Rcpp::NumericMatrix mdraw(nd,n);
   Rcpp::NumericMatrix sdraw(nd,n);
   Rcpp::IntegerVector istar(nd);
   Rcpp::NumericVector sigmadr(nd+burn);

#else

#define TRDRAW(a, b) trdraw[a][b]
#define MDRAW(a, b) mdraw[a][b]
#define SDRAW(a, b) sdraw[a][b]

void cdpmwbart(
   size_t n,
   size_t p,
   size_t np,
   double* x,
   double* y,
   double* xp,
   size_t nd,
   size_t burn,
   size_t m,
   size_t nc,
   double mybeta,
   double alpha,
   double tau,
   double mstart,
   double sstart,
   std::vector<double> &gm,
   double am,
   std::vector<double> &pm,
   std::vector<double> &gs,
   double as,
   std::vector<double> &ps,
   double lambda,
   double nu,
   size_t drawmu,
   size_t drawsd,
   double* iw,
   double iq,
   unsigned int n1, // additional parameters needed to call from C++
   unsigned int n2,
   double* sigmadr,
   double* _trdraw,
   double* trmean,
   double* _mdraw,
   double* _sdraw,
   int* istar
)
{

   //random number generation
   arn gen(n1, n2);

   size_t mgsize = gm.size();
   size_t sgsize = gs.size();

   //size_t nkeeptrain=nd;
   //std::vector<double*> mdraw(nkeeptrain), sdraw(nkeeptrain), trdraw(nkeeptrain);
   std::vector<double*> mdraw(nd), sdraw(nd), trdraw(nd);

   for(size_t i=0; i<nd; ++i) {
     size_t j=i*n;
     mdraw[i]=&_mdraw[j];
     sdraw[i]=&_sdraw[j];
     trdraw[i]=&_trdraw[j];
   }

#endif

   double ip=1.-iq;

   size_t onesigma = 1-drawsd;

   for(int i=0;i<n;i++) trmean[i]=0.0;

   //--------------------------------------------------
   printf("*****Into main of DPMBART\n");
   //--------------------------------------------------
   //check args
   printf("****p,n: %ld, %ld\n",p,n);
   printf("****np: %ld\n",np);
   printf("****nd,burn,ntree (m), nc: %ld, %ld, %ld, %ld\n",nd,burn,m,nc);
   printf("****y1,n: %lf %lf\n",y[0],y[n-1]);
   printf("****power, base, tau: %lf, %lf, %lf\n",mybeta,alpha,tau);
   printf("****mu starting value, sigma starting value: %lf, %lf\n",mstart,sstart);

   printf("****grid size and alpha for mu: %ld, %lf\n",mgsize,am);
   printf("****mu grid,  first and last: %lf, %lf\n",gm[0],gm[mgsize-1]);
   printf("****mu prior, first and last: %lf, %lf\n",pm[0],pm[mgsize-1]);

   printf("****grid size and alpha for sigma: %ld, %lf\n",sgsize,as);
   printf("****sigma grid,  first and last: %lf, %lf\n",gs[0],gs[sgsize-1]);
   printf("****sigma prior, first and last: %lf, %lf\n",ps[0],ps[sgsize-1]);
   printf("****lambda, nu: %lf, %lf\n",lambda,nu);

   printf("****mstart,sstart: %lf, %lf\n",mstart,sstart);

   printf("****onesigma: %ld\n",onesigma);
   printf("****drawmu: %ld\n",drawmu);

   //--------------------------------------------------
   //sigma,mu: places to store draws of mu_i and sigma_i
   double *svec = new double[n];
   double *mvec = new double[n];
   for(size_t i=0;i<n;i++) {
     svec[i]=iw[i]*sstart;
     mvec[i]=mstart;
   }

   //--------------------------------------------------
   heterbart bm(m);
   bm.setprior(alpha,mybeta,tau);
   double *ytemp = new double[n]; //will be y_i - mu_i on each draw
   bm.setdata(p,n,x,ytemp,nc);

   //--------------------------------------------------
   dpm mu(n,mstart);
   double *ymu = new double[n];
   mu.set_data(ymu);
   mu.set_prior(mgsize,&gm[0],&pm[0]);
   mu.set_alpha(am);
   //mu.set_sigma(1.0); //ymu will be standardized to make this right
   if(onesigma) {
      //mu.set_etaconst(sstart);
      mu.set_eta(svec);
   } else {
      mu.set_eta(svec);
   }
   mu.toscreen();

   //--------------------------------------------------
   dps sigma(n,sstart);
   double *ysigma = new double[n];
   sigma.set_data(ysigma);
   sigma.set_prior(sgsize,&gs[0],&ps[0]);
   sigma.set_alpha(as);
   sigma.toscreen();

   gen.set_df(n+nu);



   //--------------------------------------------------
   //mcmc
   printf("\nMCMC\n");
   size_t printevery=100;

   printf("\nMCMC All\n");
   time_t tp;
   int time1 = time(&tp);
   double rss=0,restemp=0;
   for(size_t i=0;i<(nd+burn);i++) {
      if(i%printevery==0) printf("done %zu (out of %lu)\n",i,nd+burn);

      //draw bart
      for(size_t k=0;k<n;k++) ytemp[k] = y[k]-mvec[k];
      bm.draw(svec,gen);

      //draw mu
      if(drawmu) {
         //if(onesigma) mu.set_etaconst(svec[0]);
         for(size_t k=0;k<n;k++) ymu[k] = (y[k]-bm.f(k));
         mu.draw(gen);
         std::vector<double> tempm = mu.theta_vector();
         for(size_t k=0;k<n;k++) mvec[k]=tempm[k];
      } else {
         for(size_t k=0;k<n;k++) mvec[k]=mstart;
      }

      //draw sigma
      if(onesigma) {
         rss=0.0;
         for(size_t k=0;k<n;k++) { 
	   restemp=(y[k] - mvec[k] - bm.f(k))/iw[k]; 
	   rss += restemp*restemp;
	   iw[k]=(restemp>0) ? 0.5/iq : 0.5/ip;
	 }
         sigmadr[i] = sqrt((nu*lambda + rss)/gen.chi_square());
         for(size_t k=0;k<n;k++) svec[k]=iw[k]*sigmadr[i];
      } else {
         for(size_t k=0;k<n;k++) {
	   ysigma[k] = (y[k] - mvec[k] - bm.f(k))/iw[k];
	   iw[k]=(ysigma[k]>0) ? 0.5/iq : 0.5/ip;
	 }
         sigma.draw(gen);
         std::vector<double> temps = sigma.theta_vector();
         for(size_t k=0;k<n;k++) svec[k]=iw[k]*temps[k];
      }

      if(i>=burn) {
	 size_t j=i-burn;
         for(size_t k=0;k<n;k++) {
	   trmean[k]+=bm.f(k);
	   TRDRAW(j,k)=bm.f(k);
	   MDRAW(j,k)=mvec[k];
	   SDRAW(j,k)=svec[k];
	 }
         istar[j]=mu.istar();
      }
   }
   int time2 = time(&tp);
   printf("time: %d\n",time2-time1);
   for(size_t k=0;k<n;k++) trmean[k]/=nd;

   if(svec) delete [] svec;
   if(mvec) delete [] mvec;
   if(ytemp) delete [] ytemp;
   if(ymu) delete [] ymu;
   if(ysigma) delete [] ysigma;

#ifndef NoRcpp
   Rcpp::List ret;
   ret["yhat.train.mean"]=trmean;
   ret["yhat.train"]=trdraw;
   ret["mdraw"]=mdraw;
   ret["sdraw"]=sdraw;
   ret["istar"]=istar;
   ret["sigmadr"]=sigmadr;

   return ret;
#endif

}
