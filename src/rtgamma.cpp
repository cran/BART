/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2019 Rodney Sparapani
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

#include "rtgamma.h"

#ifndef NoRcpp

RcppExport SEXP crtgamma(SEXP n, SEXP shape, SEXP rate, SEXP a) {
  arn gen;
  size_t N = Rcpp::as<int>(n);
//double arg1=Rcpp::as<double>(shape), arg2=Rcpp::as<double>(rate), 
//  arg3=Rcpp::as<double>(a);
  Rcpp::NumericVector z(N), A(shape), B(rate), C(a);
  size_t nA=A.size(), nB=B.size(), nC=C.size();
  for(size_t i=0; i<N; ++i) z[i]=rtgamma(A[i%nA], B[i%nB], C[i%nC], gen);
//for(size_t i=0; i<N; ++i) z[i]=rtgamma(arg1, arg2, arg3, gen);
  return Rcpp::wrap(z);
}

#endif

double rtgamma(double shape, double rate, double a, rn& gen) {
// left truncated gamma: shape must be greater than 1!
// Random Number Generation and Monte Carlo Methods, Second Edition, pp. 180-1
    if(shape<=1.) return -1.;
    double y, x=0., c=1., a_scale=a*rate, shape_shift=shape-1.,
    lambda=0.5*(a_scale-shape+sqrt(pow(a_scale-shape, 2.)+4.*a_scale))/a_scale,
      lambda_shift=1.-lambda, 
      C=1.+log(lambda_shift/shape_shift);
    while(c>x) { // do at least once
      x=gen.exp();
      y=a_scale+gen.exp()/lambda; 
      c=lambda_shift*y-shape_shift*(log(y)+C);
    }
    return y/rate;
  } 
