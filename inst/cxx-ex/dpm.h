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

#ifndef GUARD_DPM
#define GUARD_DPM

#include "dp.h"

//y_i ~ N(y; theta_i = mu, sigma_i^2), that is, DP on the mu. sigma = eta.
class dpm: public DP {
public:
   //--------------------------------------------------
   //const/dest
   dpm(): DP() {}
   dpm(size_t _n, double _theta): DP(_n,_theta) {}
   ~dpm() {}

   //--------------------------------------------------
   //public
   void toscreen() {cout << "###dpm object:\n"; DP::toscreen();}
   //the virtuals
   double f(double y, double theta, double eta);
   double draw_one_theta(std::list<size_t>& ind, rn& gen);

   //--------------------------------------------------
   //private
};


#endif
