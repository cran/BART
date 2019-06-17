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

#ifndef GUARD_common_h
#define GUARD_common_h

#ifdef MATHLIB_STANDALONE
#define NoRcpp
#else
#define RNG_Rcpp
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using std::endl;

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef NoRcpp

#include <stdio.h> // for printf

using std::cout;

#define PI 3.141592653589793238462643383280

#else // YesRcpp

#include <Rcpp.h>

#define printf Rprintf
#define cout Rcpp::Rcout

#endif

// log(2*pi)
#define LTPI 1.837877066409345483560659472811
// sqrt(2*pi)
#define RTPI 2.506628274631000502415765284811

#include "rn.h"

#endif
