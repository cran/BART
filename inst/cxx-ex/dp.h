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

#ifndef GUARD_DP
#define GUARD_DP

#include <cmath>
#include <list>
#include <algorithm>

#include "common.h"

#define RTPI 2.506628

//--------------------------------------------------------------------------------
// "theta group" class, has a theta value and indices of which obs use this theta
struct theta_gp {
public:
   //--------------------------------------------------
   //typedef
   //index iterator
   typedef std::list<size_t>::iterator iiter;

   //--------------------------------------------------
   //constructors/dest
   theta_gp():theta(0.0),ind() {}
   theta_gp(double _theta):theta(_theta),ind() {}
   ~theta_gp() {}

   //--------------------------------------------------
   // public functions
   void toscreen();

   //--------------------------------------------------
   //data members
   double theta; //the "theta" value
   std::list<size_t> ind;  //which observations have this theta value
};

//--------------------------------------------------------------------------------
// the dp class, parameter is a list of theta_gp's
class DP {
public:
   //----------------------------------------
   //typedef
   //iterate over theta clusers (theta_gp instances)
   typedef std::list<theta_gp>::iterator titer;

   //----------------------------------------
   //constructors/dest
   DP();
   DP(size_t _n, double _theta);
   virtual ~DP() {}

   //----------------------------------------
   //get,set
   void set_alpha(double alpha) {this->alpha = alpha;}
   double get_alpha()  {return alpha;}
   void set_data(double *y) {this->y = y;} //check length of y not clash with n?
   void set_prior(size_t m, double *gd, double *pri)
      {this->m=m;this->gd=gd;this->pri=pri;}

   void set_eta(double *eta) {this->eta = eta;} //check length of y not clash with n?
   double get_etaconst() {return etaconst;}
   void set_etaconst(double etaval) {etaconst =etaval;}

   //----------------------------------------
   //public functions
   size_t istar()   {return theta.size();}
   std::vector<double> theta_vector();
   std::vector<double> thetast_vector();
   void set(size_t vind, double v);
   titer find_theta(double val);
   std::vector<size_t> counts();

   size_t drind(double sum, std::vector<double>& pv, rn& gen);

   bool check();
   virtual void toscreen();

   //--------------------------------------------------
   // the two key model defining functions f and draw_one_theta, these are for prior
   // put in ones for the prior, instead of using pure virtuals
   virtual double f(double y, double theta, double eta);
   virtual double draw_one_theta(std::list<size_t>& ind, rn& gen);

   //--------------------------------------------------
   // draws
   double qo(double y, double eta);
   void draw_theta(rn& gen);
   void remix(rn& gen);
   void draw(rn& gen) { remix(gen); draw_theta(gen);}

protected:
   //--------------------
   // DP parameter
   double alpha; 		//precision (or total mass) parameter
   //note the base measure will be defined in the virtual functions
   //--------------------
   //data
   size_t n;
   double *y;                  //observations
   //--------------------
   // eta, the other parameter taken as given f(y|theta,eta), should point to of length n
   double *eta;
   double etaconst;  //if eta not set just plug in etaconst
   //--------------------
   //theta grid and prior
   size_t m;                   //size of theta grid
   double *gd;                 //theta grid
   double *pri;                //prior vals on grid
   //parameter
   std::list<theta_gp> theta;
};

void summary(DP& dpo);
void summary_cut(DP dpo);
void summary_p(DP *dpo);
#endif
