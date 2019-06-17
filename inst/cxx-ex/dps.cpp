#include "dps.h"
//################################################################################
//### f(y,theta)
double dps::f(double y, double theta, double eta) 
{
   //ignore eta, mu=0
   double r = y/theta;
   return exp(-.5*r*r)/(theta*RTPI);
}
//################################################################################
//###  dps::draw_one_theta
double dps::draw_one_theta(std::list<size_t>& ind, rn& gen) 
{
   size_t nn = ind.size();
   double s=0.0;
   double yy;
   for(theta_gp::iiter i=ind.begin();i!=ind.end();i++) {
      yy = y[*i];
      s += yy*yy;
   }

   std::vector<double> pv(m,0.0);

   for(size_t i=0;i<m;i++) {
      pv[i] = -((1.0*nn)*log(gd[i])) -s/(2.0*gd[i]*gd[i]) + log(pri[i]);
   }

   double maxel = pv[0];
   for(size_t i=1;i<m;i++) if(pv[i]>maxel) maxel = pv[i];
   double sum = 0.0;
   for(size_t i=0;i<m;i++) {
      pv[i] = exp(pv[i]-maxel);
      sum += pv[i];
   }

   size_t ii=drind(sum,pv,gen);
   return gd[ii];
}

