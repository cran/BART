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

#include "heterbartfuns.h"

//--------------------------------------------------
//heterlh, replacement for lil that only depends on sum y.
double heterlh(double b, double M, double tau) {
   double t2 =tau*tau;
   double k = b*t2+1;
   return -.5*log(k)+.5*M*M*t2/k;
}
//--------------------------------------------------
//compute b and M  for left and right give bot and v,c
void hetergetsuff(tree& x, tree::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl, double& bl, double& Ml, size_t& nr,  double& br, double& Mr, double *sigma)
{
   double *xx;//current x
   bl=0; Ml=0.0; br=0; Mr=0.0; nl=0; nr=0;
   double w;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node
         w= 1.0/(sigma[i]*sigma[i]);
         if(xx[v] < xi[v][c]) {
               nl+=1;
               bl+=w;
               Ml += w*di.y[i];
          } else {
               nr+=1;
               br+=w;;
               Mr += w*di.y[i];
          }
      }
   }
}
//--------------------------------------------------
//compute b and M for left and right bots
void hetergetsuff(tree& x, tree::tree_p l, tree::tree_p r, xinfo& xi, dinfo& di, double& bl, double& Ml, double& br, double& Mr, double *sigma)
{

   double *xx;//current x
   bl=0; Ml=0.0; br=0; Mr=0.0;
   double w;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      tree::tree_cp bn = x.bn(xx,xi);
      if(bn==l) {
         w = 1.0/(sigma[i]*sigma[i]);
         bl+=w;
         Ml += w*di.y[i];
      }
      if(bn==r) {
         w = 1.0/(sigma[i]*sigma[i]);
         br+=w;
         Mr += w*di.y[i];
      }
   }
}
//--------------------------------------------------
//draw one mu from post
double heterdrawnodemu(double b, double M, double tau, rn& gen)
{
   double muhat = M/b;
   double a = 1.0/(tau*tau);
   return (b*muhat)/(a+b) + gen.normal()/sqrt(a+b);
}
//--------------------------------------------------
//get sufficients stats for all bottom nodes, this way just loop through all the data once.
void heterallsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<double>& bv, std::vector<double>& Mv,double *sigma)
{
   tree::tree_cp tbn; //the pointer to the bottom node for the current observations
   size_t ni;         //the  index into vector of the current bottom node
   double *xx;        //current x

   bnv.clear();
   x.getbots(bnv);

   typedef tree::npv::size_type bvsz;
   bvsz nb = bnv.size();
   bv.resize(nb);
   Mv.resize(nb);

   std::map<tree::tree_cp,size_t> bnmap;
   for(bvsz i=0;i!=bnv.size();i++) {bnmap[bnv[i]]=i;bv[i]=0;Mv[i]=0.0;}

   double w;
   for(size_t i=0;i<di.n;i++) {
      w = 1.0/(sigma[i]*sigma[i]);
      xx = di.x + i*di.p;
      tbn = x.bn(xx,xi);
      ni = bnmap[tbn];

      bv[ni] += w;
      Mv[ni] += w*di.y[i];
   }
}
//--------------------------------------------------
//heter version of drmu, need b and M instead of n and sy
void heterdrmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, rn& gen)
{
   tree::npv bnv;
   std::vector<double> bv;
   std::vector<double> Mv;
   heterallsuff(t,xi,di,bnv,bv,Mv,sigma);
   for(tree::npv::size_type i=0;i!=bnv.size();i++)
      bnv[i]->settheta(heterdrawnodemu(bv[i],Mv[i],pi.tau,gen));
}
