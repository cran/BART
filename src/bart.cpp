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

#include "bart.h"

//--------------------------------------------------
//constructor
bart::bart():m(200),t(m),pi(),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),di(),dartOn(false) {}
bart::bart(size_t im):m(im),t(m),pi(),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),di(),dartOn(false) {}
bart::bart(const bart& ib):m(ib.m),t(m),pi(ib.pi),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),di(),dartOn(false)
{
   this->t = ib.t;
}
bart::~bart()
{
   if(allfit) delete[] allfit;
   if(r) delete[] r;
   if(ftemp) delete[] ftemp;
}

//--------------------------------------------------
//operators
bart& bart::operator=(const bart& rhs)
{
   if(&rhs != this) {

      this->t = rhs.t;
      this->m = t.size();

      this->pi = rhs.pi;

      p=0;n=0;x=0;y=0;
      xi.clear();

      if(allfit) {delete[] allfit; allfit=0;}
      if(r) {delete[] r; r=0;}
      if(ftemp) {delete[] ftemp; ftemp=0;}

   }
   return *this;
}
//--------------------------------------------------
//get,set
void bart::setm(size_t m)
{
   t.resize(m);
   this->m = t.size();

   if(allfit && (xi.size()==p)) predict(p,n,x,allfit);
}

//--------------------------------------------------
void bart::setxinfo(xinfo& _xi)
{
   size_t p=_xi.size();
   xi.resize(p);
   for(size_t i=0;i<p;i++) {
     size_t nc=_xi[i].size();
      xi[i].resize(nc);
      for(size_t j=0;j<nc;j++) xi[i][j] = _xi[i][j];
   }
}
//--------------------------------------------------
void bart::setdata(size_t p, size_t n, double *x, double *y, size_t numcut)
{
  int* nc = new int[p];
  for(size_t i=0; i<p; ++i) nc[i]=numcut;
  this->setdata(p, n, x, y, nc);
  delete [] nc;
}

void bart::setdata(size_t p, size_t n, double *x, double *y, int *nc)
{
   this->p=p; this->n=n; this->x=x; this->y=y;
   if(xi.size()==0) makexinfo(p,n,&x[0],xi,nc);

   if(allfit) delete[] allfit;
   allfit = new double[n];
   predict(p,n,x,allfit);

   if(r) delete[] r;
   r = new double[n];

   if(ftemp) delete[] ftemp;
   ftemp = new double[n];

   di.n=n; di.p=p; di.x = &x[0]; di.y=r;
   for(size_t j=0;j<p;j++){
     nv.push_back(0);
     pv.push_back(1/(double)p);
   }
}
//--------------------------------------------------
void bart::predict(size_t p, size_t n, double *x, double *fp)
//uses: m,t,xi
{
   double *fptemp = new double[n];

   for(size_t j=0;j<n;j++) fp[j]=0.0;
   for(size_t j=0;j<m;j++) {
      fit(t[j],xi,p,n,x,fptemp);
      for(size_t k=0;k<n;k++) fp[k] += fptemp[k];
   }

   delete[] fptemp;
}
//--------------------------------------------------
void bart::draw(double sigma, rn& gen)
{
   for(size_t j=0;j<m;j++) {
      fit(t[j],xi,p,n,x,ftemp);
      for(size_t k=0;k<n;k++) {
         allfit[k] = allfit[k]-ftemp[k];
         r[k] = y[k]-allfit[k];
      }
      bd(t[j],xi,di,pi,sigma,nv,pv,aug,gen);
      drmu(t[j],xi,di,pi,sigma,gen);
      fit(t[j],xi,p,n,x,ftemp);
      for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
   }
   if(dartOn) draw_s(gen);
}
//--------------------------------------------------
void bart::draw_s(rn& gen)
{
  // Draw sparsity parameter theta (Linero calls it alpha); see Linero, 2018
  // theta / (theta + rho ) ~ Beta(a,b)
  // Set (a=0.5, b=1) for sparsity
  // Set (a=1, b=1) for non-sparsity
  // rho = p usually, but making rho < p increases sparsity
  if(!const_theta){
  double sumlogpv=0.,mx,sm,log_sum_exp;
  std::vector<double> lambda_g (1000,0.);  
  std::vector<double> theta_g (1000,0.);
  std::vector<double> lwt_g (1000,0.);
    for(size_t j=0;j<p;j++) sumlogpv+=log(pv[j]);
    for(size_t k=0;k<1000;k++){
      lambda_g[k]=(double)(k+1)/1001.;
      theta_g[k]=(lambda_g[k]*rho)/(1.-lambda_g[k]);
      lwt_g[k]=lgamma(theta_g[k])-lgamma(theta_g[k]/(double)p)*(double)p
	       +(theta_g[k]/(double)p)*sumlogpv
	       +(a-1)*log(lambda_g[k])+(b-1)*log(1-lambda_g[k]);
      if(k==0) mx=lwt_g[0];
      else if(lwt_g[k]>mx) mx=lwt_g[k];
    }
    //for(size_t k=0;k<1000;k++) if(lwt_g[k]>mx) mx=lwt_g[k];
    sm=0.;
    for(size_t k=0;k<1000;k++) sm+=exp(lwt_g[k]-mx);
    log_sum_exp=mx+log(sm);
    for(size_t k=0;k<1000;k++) lwt_g[k]=exp(lwt_g[k]-log_sum_exp);
    gen.set_wts(lwt_g);
    theta=theta_g[gen.discrete()];
   } 
  // Now draw s, the vector of splitting probabilities
  std::vector<double> _alpha(p);
  for(size_t j=0;j<p;j++) _alpha[j]=theta/(double)p+nv[j];
  gen.set_alpha(_alpha);
  pv=gen.dirichlet();
}
//--------------------------------------------------
//public functions
void bart::pr() //print to screen
{
   cout << "*****bart object:\n";
   cout << "m: " << m << std::endl;
   cout << "t[0]:\n " << t[0] << std::endl;
   cout << "t[m-1]:\n " << t[m-1] << std::endl;
   cout << "prior and mcmc info:\n";
   pi.pr();
   if(dart){
     cout << "*****dart prior (On):\n";
     cout << "a: " << a << std::endl;
     cout << "b: " << b << std::endl;
     cout << "rho: " << rho << std::endl;
     cout << "augmentation: " << aug << std::endl;
   }
   else cout << "*****dart prior (Off):\n";
   if(p) cout << "data set: n,p: " << n << ", " << p << std::endl;
   else cout << "data not set\n";
}
