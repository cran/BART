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

#include "dp.h"

//###--------------------------------------------------
//### DP constructors
DP::DP():alpha(1.0),n(0),y(nullptr),eta(nullptr),etaconst(0),m(0),gd(nullptr),pri(nullptr),theta() {}
//--------------------------------------------------
//inititialize to one group with the input theta value
DP::DP(size_t _n, double _theta) : alpha(1.0),n(0),y(nullptr),eta(nullptr),etaconst(0),m(0),gd(nullptr),pri(nullptr),theta()
{
   n = _n;

   //make a single theta group
   theta_gp th;
   th.theta = _theta;
   th.ind.clear();

   //assign all observations to theta group
   for(size_t i=0;i<n;i++) {
      th.ind.push_back(i);
   }

   //add the single group to the theta list
   theta.push_back(th);
}
//###--------------------------------------------------
//### DP to screen
void DP::toscreen()
{
   cout << "*****DP object:\n";
   cout << "\t alpha: " << alpha << std::endl;
   cout << "\t n: " << n << std::endl;
   cout << "\t number of clusters: " << theta.size() << endl;
   for(titer i=theta.begin();i!=theta.end();i++)
      cout << "\t\t theta: " << i->theta <<  "; cluster size:" << i->ind.size() << endl;
   cout << "\t m: (size of grid): " << m << endl;
   if(y==nullptr) {
      cout << "\t data not set\n";
   } else {
      cout << "\t y: " << y[0] << ", " << y[n-1] << endl;
   }
   if(eta==nullptr) {
      cout << "\t eta not set\n";
   } else {
      cout << "\t eta: " << eta[0] << ", " << eta[n-1]  << endl;
   }
   cout << "etaconst: " << etaconst << endl;
   cout << "\n\n";
}
//###--------------------------------------------------
//### get theta vector with value for each observation (so some may be repeated)
std::vector<double> DP::theta_vector()
{
   std::vector<double> v(n);
   for(titer i=theta.begin(); i!=theta.end();i++) {
      for(theta_gp::iiter j=i->ind.begin();j!=i->ind.end();j++) v[*j] = i->theta;
   }
   return v;
}
//###--------------------------------------------------
//### get theta star vector
std::vector<double> DP::thetast_vector()
{
   std::vector<double> tst;
   for(titer i= theta.begin();i!=theta.end();i++) tst.push_back(i->theta);
   return tst;
}

//###--------------------------------------------------
//###check that the thetastar are unique, and the the union of indices is 0,1,2,,n-1
bool DP::check()
{
   std::vector<double> tst = thetast_vector();
   std::sort(tst.begin(),tst.end());
   for(std::vector<double>::size_type i=1;i<tst.size();i++) if(tst[i]==tst[i-1]) return false;

   std::list<size_t> allind;
   for(titer i=theta.begin();i!=theta.end();i++) copy(i->ind.begin(),i->ind.end(),back_inserter(allind));
   allind.sort();

   size_t cnt=0;
   for(theta_gp::iiter i=allind.begin();i!=allind.end();i++) {
      if(*i != cnt) return false;
      //std::cout << *i << std::endl;
      cnt++;
   }
   return true;
}
//###--------------------------------------------------
//### find iterator of cluster with theta = val
DP::titer DP::find_theta(double val)
{
   titer i = theta.begin();
   bool not_found = true;
   while(not_found && i!=theta.end()) {
      if(i->theta == val) not_found = false;
      else i++;
   }
   return i;
}
//###--------------------------------------------------
//### set the theta value for a particular observation, HARD!!
// ### set theta[vind] to v
void DP::set(size_t vind, double v)
{
   //should check that n>0

   titer i  = theta.begin(); //which cluster is vind in
   theta_gp::iiter j; //which index in the cluster is vind

   bool not_found = true;
   while(not_found && i != theta.end()) {
       j = find(i->ind.begin(),i->ind.end(),vind);
       if(j != i->ind.end())
          not_found = false;
        else
          i++;
   }
   //if(not_found) {cout << "error: index not found in DP::set\n"; exit(EXIT_FAILURE);}

   if(i->theta != v) { //if it already has the value, nothing to do
      i->ind.erase(j); // drop it from old group
      if(i->ind.size() == 0) theta.erase(i); //if group now empty, drop it from theta
      titer k = find_theta(v); //see if you can find a theta with value v
      if(k == theta.end()) { //if you can't make a new group and add it to theta
         theta_gp tg;
         tg.theta = v;
         tg.ind.push_back(vind);
         theta.push_back(tg);
      } else { //if you can, just add the index to the list for the group
         k->ind.push_back(vind);
      }
   }
}
//###--------------------------------------------------
//### get the size of each theta cluster
std::vector<size_t> DP::counts()
{
   std::vector<size_t> cnts;
   for(titer i= theta.begin();i!=theta.end();i++) cnts.push_back(i->ind.size());
   return cnts;
}
//###--------------------------------------------------
//### draw from a discrete, note that sum is sum of pv, so pv may be unnormalized
size_t DP::drind(double sum,std::vector<double>& pv, rn& gen)
{
   double u = gen.uniform();
   double psum=pv[0]/sum;
   size_t ii=0;
   while(psum<u) {
      ii++;
      psum += pv[ii]/sum;
   }
   return ii;
}
//###--------------------------------------------------
//### this is just to play around with ideas in Chapter 13 of "Accelerated .."
void summary(DP& dpo) {
   dpo.toscreen();
}
//###--------------------------------------------------
//### this is just to play around with ideas in Chapter 13 of "Accelerated .."
void summary_cut(DP dpo) {
   dpo.toscreen();
}
//###--------------------------------------------------
//### this is just to play around with ideas in Chapter 13 of "Accelerated .."
void summary_p(DP *dpo) {
   dpo->toscreen();
}
//###--------------------------------------------------
//### f(y,theta)=1
double DP::f(double y, double theta, double eta)
{
   return 1.0;
}
//###--------------------------------------------------------------------------------
//### DP::draw_one_theta
double DP::draw_one_theta(std::list<size_t>& ind, rn& gen)
{

   std::vector<double> pv(m);
   for(size_t i=0;i<m;i++) pv[i] = pri[i];
   size_t ii=drind(1.0,pv,gen);
   return gd[ii];
}
//###--------------------------------------------------------------------------------
//### DP::qo(double y) the base measure prior predictive at y
double DP::qo(double y, double eta)
{
   //if(m==0) {cout << "error: index not found in DP::set\n"; exit(EXIT_FAILURE);}

   double q = 0.0;
   for(size_t i=0;i<m;i++) q += f(y,gd[i],eta)*pri[i];
   return q;
}
//###--------------------------------------------------------------------------------
//### DP::draw_theta, this does the basic dp algorithm
void DP::draw_theta(rn& gen) {

   //get counts and thetastar for initial partition
   std::vector<double> tst = thetast_vector();
   std::vector<size_t> cnt = counts();
   std::vector<size_t> icnt = cnt; //need to know the initial number
   size_t np = cnt.size(); //initial number of partitions or clusters

   //this is tricky code because dimension changes, new partitions get put on the end
   //   after you are done you will drop any empties
   //   i loops over partitions but par manages the one you are on
   //   j loops over observations but obs manages the one you are on.
   titer par = theta.begin();
   for(size_t i=0;i<np;i++) {
      theta_gp::iiter ob = par->ind.begin();
      for(size_t j=0;j<icnt[i];j++) {
         //std::cout << "doing obs: " << *ob << ", " << y[*ob] << std::endl;
	 size_t nts = tst.size();
         std::vector<double> pv(nts + 1);
	 double sum=0.0;
	 for(size_t k=0;k<nts;k++) {
	    if(cnt[k] == 0) {
               pv[k] = 0.0;
            } else {
               if(eta == nullptr) {
	          pv[k] = cnt[k]*f(y[*ob],tst[k],etaconst);
               } else {
	          pv[k] = cnt[k]*f(y[*ob],tst[k],eta[*ob]);
               }
            }
	    sum += pv[k];
	 }
         if(eta == nullptr) {
	    pv[nts] = alpha * qo(y[*ob],etaconst);
         } else {
	    pv[nts] = alpha * qo(y[*ob],eta[*ob]);
         }
	 sum += pv[nts];
	 size_t ii = drind(sum,pv,gen);

	 //for(size_t jj = 0;jj<nts;jj++) std::cout << "i,j,jj,tst[jj],pv[jj]: "
	 //              << i << " " << j << " " << jj << " " << tst[jj] << " " << (pv[jj]/sum) << std::endl;
         //std::cout << (pv[nts]/sum) << std::endl;
	 //std::cout << "ii: " << ii << std::endl;

	 if(ii==nts) {  // birth
	    std::list<size_t> tind; tind.push_back(*ob);
	    double ntheta = draw_one_theta(tind,gen); //draw new theta
	    tst.push_back(ntheta); cnt.push_back(1); cnt[i] -=1;
	    theta_gp ngp; ngp.theta = ntheta; ngp.ind.push_back(*ob);
	    theta.push_back(ngp);
	    ob = par->ind.erase(ob);
	 } else if(tst[ii] == par->theta) { //drew yourself
	    ob++;
	 } else {  //drew one of the other old thetas
	    cnt[i] -=1; cnt[ii] += 1;
	    titer tp = find_theta(tst[ii]); tp->ind.push_back(*ob);
	    ob = par->ind.erase(ob);
	 }
      }
      par++;
   }
   //drop all empty partitions
   par = theta.begin();
   while(par!=theta.end()) {
      if(par->ind.size() == 0) par = theta.erase(par);
      else  par++;
   }
}
//###--------------------------------------------------
//### DP::remix, redraw theta in each cluster.
void DP::remix(rn& gen)
{
   for(titer i = theta.begin(); i!= theta.end(); i++) {
      i->theta = draw_one_theta(i->ind,gen);
   }
}
//--------------------------------------------------
//--------------------------------------------------
//theta_gp
//--------------------------------------------------
//theta_gp to screen
void theta_gp::toscreen()
{
   cout << "theta group:\n";
   cout << "theta: " << theta << endl;
   cout << "size of partition: " << ind.size() << endl;
   //for(iiter i = ind.begin(); i!=ind.end();i++) std::cout << *i << " ";
   cout << endl;
}

