
library(BART)

## replication script for JSS article

## directory to create graphics files in, if any
## if option not specified, then default assumed to be 'NONE'
options(figures='.')
## options(figures='../vignettes/figures')

## for single-threading, specify one core
## Windows lacks forking (and generally lacks OpenMP)
## so single-threading only
## for multi-threading, specify the number of cores
## a Unix-like OS provides forking for multi-threading
## (and often OpenMP is available as well)

if(.Platform$OS.type=='unix') {
    ## there are diminishing returns so often 8 cores is sufficient
    options(mc.cores=min(8, parallel::detectCores()))
} else {
    options(mc.cores=1)
}

## uncomment these options to compare multiple threading
## with a single thread as a double-check of the package
## due to random seed/stream progression the results will
## NOT be identical, but they should be very comparable
## options(mc.cores=1, figures='../vignettes/single')

## Section 3, The Boston Housing Data including Figures 1-6
source(system.file('demo/boston.R', package='BART'))

## Section 4.2, Probit BART Example: Chronic Pain and Obesity
## Figure 7
source(system.file('demo/nhanes.pbart1.R', package='BART'))
## Figure 8
source(system.file('demo/nhanes.pbart2.R', package='BART'))

## Section 4.4, Multinomial BART Example: Alligator Food Preference
## Figure 9
source(system.file('demo/alligator.R', package='BART'))

## Section 4.5, Convergence Diagnostics for Binary and Categorical Outcomes
## Figures 10-12
source(system.file('demo/geweke.pbart2.R', package='BART'))

## Section 4.6, BART and Variable Selection
## Figure 13
source(system.file('demo/sparse.pbart.R', package='BART'))

## Section 5.1, Survival Analysis with BART Example: Advanced Lung Cancer
## Figure 14
source(system.file('demo/lung.surv.bart.R', package='BART'))

## Section 5.3, Competing Risks with BART Example: Liver Transplants
## Figure 15
source(system.file('demo/liver.crisk.bart.R', package='BART'))

## Section 5.4, Recurrent Events with BART Example: Bladder Tumors}
## Figures 16-18
source(system.file('demo/bladder.recur.bart.R', package='BART'))

figures = getOption('figures', default='NONE')

## Figure 19
library(Rcpp)

dbetapr = function (x, shape1, shape2, scale = 1, log = FALSE)
    cpp_dbetapr(x, shape1, shape2, scale, log[1L])

sourceCpp(code='
#include <Rcpp.h>
#define GETV(x, i)      x[i % x.length()]    // wrapped indexing of vector

// [[Rcpp::plugins(cpp11)]]

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using Rcpp::NumericVector;
using std::log1p;
/*
*  Beta prime distribution
*  Values:
*  x > 0
*  Parameters:
*  alpha > 0
*  beta > 0
*  sigma > 0
*/
inline double logpdf_betapr(double x, double alpha, double beta,
                            double sigma, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(alpha) || ISNAN(beta) || ISNAN(sigma))
    return x+alpha+beta+sigma;
#endif
  if (alpha <= 0.0 || beta <= 0.0 || sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x <= 0.0 || !R_FINITE(x))
    return R_NegInf;
  double z = x / sigma;
  // pow(z, alpha-1.0) * pow(z+1.0, -alpha-beta) / R::beta(alpha, beta) / sigma;
  return log(z) * (alpha-1.0) + log1p(z) * (-alpha-beta) -
    R::lbeta(alpha, beta) - log(sigma);
}

// [[Rcpp::export]]
NumericVector cpp_dbetapr(
    const NumericVector& x,
    const NumericVector& alpha,
    const NumericVector& beta,
    const NumericVector& sigma,
    const bool& log_prob = false
  ) {

  if (std::min({x.length(), alpha.length(),
                beta.length(), sigma.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    alpha.length(),
    beta.length(),
    sigma.length()
  });
  NumericVector p(Nmax);

  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_betapr(GETV(x, i), GETV(alpha, i),
                         GETV(beta, i), GETV(sigma, i),
                         throw_warning);

  if (!log_prob)
    p = Rcpp::exp(p);

  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}
')

x=seq(0, 5, length.out=1001)

plot(x, dbetapr(x, 0.5, 1, 0.5), col='gray', type='l', lty=4, lwd=2,
     log='y', xlab=expression(italic(x)),
     ylab=expression(italic(log(f(x, a, b, rho/P)))))
lines(x, dbetapr(x, 0.5, 1, 1), col='blue', lty=3, lwd=2)
lines(x, dbetapr(x, 1, 1, 0.5), col='red', lty=2, lwd=2)
lines(x, dbetapr(x, 1, 1, 1), lwd=2)
legend('topright', col=c('black', 'red', 'blue', 'gray'), lwd=2, lty=1:4,
       legend=c(expression(italic(log(f(x, 1, 1, 1)))),
                expression(italic(log(f(x, 1, 1, 0.5)))),
                expression(italic(log(f(x, 0.5, 1, 1)))),
                expression(italic(log(f(x, 0.5, 1, 0.5))))))
if(figures!='NONE')
    dev.copy2pdf(file=paste(figures, 'sparse-beta-prime.pdf', sep='/'))

##Figure 20
plot(0, 1, type='h', lwd=3,
     xlim=c(-0.5, 5.5), ylim=c(0, 1),
     ylab='Proportionate length of chain processing time', xlab='Chains')
for(i in 1:5) lines(i, 0.28, type='h', lwd=3, col=2)
lines(c(-0.25, 5.25), c(0.1, 0.1), type='l', lty=2)
text(-0.35, 0.1, labels='b')
if(figures!='NONE')
    dev.copy2pdf(file=paste(figures, 'parallel.pdf', sep='/'))

##Figure 21
C <- 1:64
for(b in c(0.025, 0.1)) {
    amdahl <- 1/(b+(1-b)/C)
    if(b==0.025)
        plot(C, amdahl, type='l', col=2,
             xlim=c(1, 80), ylim=c(0, 30), log='x',
             xlab='B: number of CPU', ylab='Gain')
    else
        lines(C, amdahl, type='l')
    text(78, amdahl[64], labels=paste0(b))
}
if(figures!='NONE')
    dev.copy2pdf(file=paste(figures, 'amdahl.pdf', sep='/'))
