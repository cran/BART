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

#include "latent.h"

#ifndef NoRcpp

RcppExport SEXP cdraw_lambda(SEXP lambda, SEXP mean, SEXP kmax, SEXP thin) {
  rk_state *state;
  newRNGstates();
  state = states[0];
  double x=draw_lambda_i(Rcpp::as<double>(lambda), Rcpp::as<double>(mean),
			 Rcpp::as<int>(kmax), Rcpp::as<int>(thin), state);
  deleteRNGstates();
  return Rcpp::wrap(x);
}

RcppExport SEXP cdraw_z(SEXP mean, SEXP tau, SEXP lambda) {
  rk_state *state;
  newRNGstates();
  state = states[0];
  double z=rtnorm_reject(Rcpp::as<double>(mean), Rcpp::as<double>(tau),
			 sqrt(Rcpp::as<double>(lambda)), state);
  deleteRNGstates();
  return Rcpp::wrap(z);
}

#endif

/*
 * draw_lambda_prior:
 *
 * draw Lambda from its (infinite mixture) prior
 */

double draw_lambda_prior(double *psii, int kmax, rk_state *state)
{
  double lambda;
  int k;

  lambda = 0.0;
  for(k=0; k<=kmax; k++) {
    lambda += psii[k] * expo_rand(state);
  }

  return lambda;
}


/*
 * draw_lambda:
 *
 * Metropolis-Hastings algorithm for drawing lambda from its
 * full conditional -- uses proposals from the prior
 */

double draw_lambda_i(double lambda_old, double xbeta,
                     int kmax, int thin, rk_state *state)
{
  int t, k;
  double lambda, lp, lpold, m, s;
  double *psii;

  /* calculate the probability og the previous lambda */
  s = sqrt(lambda_old);
  m = xbeta;
  lpold = pnorm(0.0, m, s, 0, 1);

  /* allocate psii */
  psii = (double*) malloc(sizeof(double) * (kmax+1));
  for(k=0; k<=kmax; k++) psii[k] =  2.0/((1.0+k)*(1.0+k));

  /* thinning is essential when kappa is large */
  for(t=0; t<thin; t++) {

    /* propose a new lambda from the prior */
    lambda = draw_lambda_prior(psii, kmax, state);

    /* calculate the probability of the propsed lambda */
    s = sqrt(lambda);
    m = xbeta;
    lp = pnorm(0.0, m, s, 0, 1);

    /* MH accept or reject */
    if(runi(state) < exp(lp - lpold)) {
      lambda_old = lambda;
      lpold = lp;
    }
  }

  /* possibly clean up psii */
  free(psii);

  return lambda_old;
}


/*
 * draw_lambda:
 *
 * R interface to lambda draw from the posterior conditional
 * by rejection sampling and proposals from the prior
 */

void draw_lambda(int n, double *xbeta_in, int kmax, 
                   int thin, double *lambda_inout)
{
#ifdef _OPENMP
  #pragma omp parallel
  {
  int i, start, step;
  rk_state *state;
  start = omp_get_thread_num();
  step = omp_get_max_threads();
#else
  int i, start, step;
  rk_state *state;
  start = 0;
  step = 1;
#endif

  state = states[start];

  for(i=start; i<n; i+=step) {
      lambda_inout[i] = draw_lambda_i(lambda_inout[i], xbeta_in[i],
                                      kmax, thin, state);
  }

#ifdef _OPENMP
  }
#endif
}


/*
 * draw_z:
 *
 * C-side for R function draw.z, which draws from the full
 * conditional of the latent z-values given y, beta, and
 * the latent lambdas
 */

void draw_z(int n, double *xbeta, double *lambda, double *z_out)
{
#ifdef _OPENMP
  #pragma omp parallel
  {
  int i, start, step;
  double aux[2];
  double lambdai_sqrt, xbetai;
  rk_state *state;
  start = omp_get_thread_num();
  step = omp_get_max_threads();
#else
  int i, start, step;
  double aux[2];
  double lambdai_sqrt, xbetai;
  rk_state *state;
  start = 0;
  step = 1;
#endif

  state = states[start];

  /* loop over rows of X */
  for(i=start; i<n; i+=step) {

    /* calculate the mean and variance of the normal */
    xbetai = xbeta[i];
    lambdai_sqrt = sqrt(lambda[i]);

    /* draw until we get one in the right half-plane */
    if(xbetai >= 0) {
      do {
        rnor(aux, state);
        z_out[i] = aux[0]*lambdai_sqrt + xbetai;
        if(z_out[i] >= 0) break;
        z_out[i] = aux[1]*lambdai_sqrt + xbetai;

      } while (z_out[i] < 0.0);
    } else {
      z_out[i] = rtnorm_reject(xbetai, 0.0, lambdai_sqrt, state);
      assert(z_out[i] > 0.0);
    }
  }
#ifdef _OPENMP
  }
#endif
}
