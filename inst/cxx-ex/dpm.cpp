#include "dpm.h"

//################################################################################
//### f(y,theta)
double dpm::f(double y, double theta, double eta) 
{
   double sigma = eta;
   double r = (y-theta)/sigma;
   return exp(-.5*r*r)/(sigma*RTPI);
}
//################################################################################
//###  dpm::draw_one_theta
double dpm::draw_one_theta(std::list<size_t>& ind, rn& gen)
{
   double mmu,mvar; //mean and var for suff stat normal likelihood

   if(eta == nullptr) {
      double sigma = etaconst;
   
      size_t nn = ind.size();
      double ybar=0.0;
      double yy;
      for(theta_gp::iiter i=ind.begin();i!=ind.end();i++) {
         yy = y[*i];
         ybar += yy;
      }
      mmu = ybar/(1.0*nn);
      mvar = (sigma*sigma)/(1.0*nn);
   } else {
      mvar = 0.0;
      mmu = 0.0;
      for(theta_gp::iiter i=ind.begin();i!=ind.end();i++) {
         double sigma = eta[*i];
         double w = 1.0/(sigma*sigma);
         mvar += w;
         mmu += w*y[*i];
      }
      mvar = 1.0/mvar;
      mmu *= mvar;
   }
   /*
   std::cout << "ybar: " << ybar << std::endl;
   std::cout << "sigma: " << sigma << std::endl;
   std::cout << "nn: " << nn << std::endl;
   */

   std::vector<double> pv(m,0.0);

   for(size_t i=0;i<m;i++) {
      //pv[i] = -(nn*log(gd[i])) -s/(2.0*gd[i]*gd[i]) + log(pri[i]);
      //double r = (ybar-gd[i])/sigma;
      //std::cout << "r: " << r << " ,, "; 
      //pv[i] = -(1.0*nn)*r*r/2.0 + log(pri[i]);
      pv[i] = -.5*(gd[i]-mmu)*(gd[i]-mmu)/mvar + log(pri[i]);
   }
   /*
   for(size_t i=0;i<m;i++) {
      std::cout << pv[i] << ", ";
   }
   */

   double maxel = pv[0];
   for(size_t i=1;i<m;i++) if(pv[i]>maxel) maxel = pv[i];
   double sum = 0.0;
   for(size_t i=0;i<m;i++) {
      pv[i] = exp(pv[i]-maxel);
      sum += pv[i];
   }
   /*
   for(size_t i=0;i<m;i++) {
      std::cout << pv[i] << ", ";
   }
   */

   size_t ii=drind(sum,pv,gen);
   return gd[ii];
}

