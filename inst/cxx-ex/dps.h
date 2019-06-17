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

#ifndef GUARD_DPS
#define GUARD_DPS

#include "dp.h"

//y_i ~ N(y; mu=0, (theta_i=sigma)^2), that is, DP on the sigma.
class dps: public DP {
public:
   //--------------------------------------------------
   //const/dest
   dps(): DP() {}
   dps(size_t _n, double _theta): DP(_n,_theta) {}
   ~dps() {}

   //--------------------------------------------------
   // let user know there is no eta.
   void set_eta(double *eta) {cout << "there is no eta for this model\n";}
   double get_etaconst() {cout << "there is no eta for this model\n"; return 0.0;}
   void set_etaconst(double eval) {cout << "there is no eta for this model\n";}

   //--------------------------------------------------
   //public
   void toscreen() {cout << "###dps object:\n"; DP::toscreen();}
   //the virtuals
   double f(double y, double theta, double eta);
   double draw_one_theta(std::list<size_t>& ind, rn& gen);
};


#endif
