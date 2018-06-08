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

#ifndef GUARD_latent
#define GUARD_latent

#include "common.h"
#include "randomkit.h"
#include <stdlib.h>
#include <assert.h>
#include "rand_draws.h"

#ifdef NoRcpp
using ::pnorm;
#else
using R::pnorm;
#endif

extern int NS;
extern rk_state** states;

#ifndef NoRcpp

RcppExport SEXP cdraw_lambda(SEXP lambda, SEXP mean, SEXP kmax, SEXP thin);

RcppExport SEXP cdraw_z(SEXP mean, SEXP tau, SEXP lambda);

#endif

void draw_z(int n, double *xbeta, double *lambda, double *z_out);
void draw_lambda(int n, double *xbeta_in, int kmax, int thin, double *lambda_inout);
double draw_lambda_i(double lambda_old, double xbeta, int kmax, int thin, rk_state *state);
double draw_lambda_prior(double *psii, int kmax, rk_state *state);

#endif
