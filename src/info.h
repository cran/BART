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

#ifndef GUARD_info_h
#define GUARD_info_h
#include "common.h"
//data
class dinfo {
public:
   dinfo() {p=0;n=0;x=0;y=0;}
   size_t p;  //number of vars
   size_t n;  //number of observations
   double *x; // jth var of ith obs is *(x + p*i+j)
   double *y; // ith y is *(y+i) or y[i]
};
//prior and mcmc
class pinfo
{
public:
   pinfo(): pbd(1.0),pb(.5),alpha(.95),mybeta(2.0),tau(1.0) {}
//mcmc info
   double pbd; //prob of birth/death
   double pb;  //prob of birth
//prior info
   double alpha;
   double mybeta;
   double tau;
   void pr() {
      cout << "pbd,pb: " << pbd << ", " << pb << std::endl;
      cout << "alpha,beta,tau: " << alpha << 
             ", " << mybeta << ", " << tau << std::endl;
   }
};

#endif
