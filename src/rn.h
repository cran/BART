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

#ifndef GUARD_rn_h
#define GUARD_rn_h

//pure virtual base class for random numbers
class rn
{
 public:
  rn() {}
  virtual double normal() = 0; //standard normal
  virtual double uniform() = 0; //uniform(0,1)
  virtual double chi_square() = 0; //chi-square
  virtual double exp() = 0; //exponential
  virtual double gamma() = 0; //gamma
  virtual size_t discrete() = 0; //discrete (categorical) distribution
  virtual size_t geometric() = 0; //geometric distribution
  virtual void set_df(double _df) = 0; //set df for chi-square
  virtual double get_df() = 0;
  virtual void set_shape(double _shape) = 0; //set shape parameter for gamma
  virtual double get_shape() = 0;
  virtual void set_wts(std::vector<double>& _wts) = 0;
  virtual double get_p() = 0;
  virtual void set_p(double _p) = 0;
  virtual double beta() = 0; 
  virtual void set_beta(double _a, double _b) = 0; 
  virtual double get_a() = 0;
  virtual double get_b() = 0;
  virtual std::vector<double> dirichlet() = 0; 
  virtual void set_alpha(std::vector<double>& _alpha) = 0;
  virtual std::vector<double> get_alpha() = 0;
  virtual ~rn() {}
};

#ifdef RNG_random

#include <Rmath.h>
#include <random> 

//abstract random number generator based on C++ <random>
class arn: public rn
{
  //typedefs
  typedef std::default_random_engine genD;
  typedef std::normal_distribution<double> norD;
  typedef std::uniform_real_distribution<double> uniD;
  typedef std::chi_squared_distribution<double> chiD;
  typedef std::gamma_distribution<double> gamD;
//  typedef std::discrete_distribution<int> disD;
  typedef std::geometric_distribution<int> geoD;
 public:
  //constructor
  arn() {}
  arn(unsigned int n1, unsigned int n2) {this->n1=n1; this->n2=n2;}
  //virtual
  virtual ~arn() {}
  virtual double normal() {return (nor)(gen);}
  virtual double uniform() {return (uni)(gen);}
  virtual double chi_square() {return (chi)(gen);}
  virtual double exp() {return -log(uniform());}
  virtual double gamma() {return (gam)(gen);}
//  virtual size_t discrete() {return (dis)(gen);}
  virtual size_t discrete() {
    size_t p=wts.size(), x=0;
    std::vector<int> vOut (p,0);
    ::rmultinom(1,&wts[0],p,&vOut[0]); 
    if(vOut[0]==0) for(size_t j=1;j<p;j++) x += j*vOut[j]; 
    return x;
  } 
  virtual size_t geometric() {return (geo)(gen);}
  virtual void set_df(double _df){
    this->df=_df; // set class-member df
    chi=chiD(_df); // assign df to chi-square dist object
  }
  virtual double get_df()  {return df;}
  virtual void set_shape(double _shape){
    this->shape=_shape; // set class-member shape
    gam=gamD(_shape,1.); // assign shape to gamma dist object
  }
  virtual double get_shape() {return shape;}
  virtual void set_wts(std::vector<double>& _wts){
    this->wts=_wts;
//    dis=disD(wts);
  }
  virtual double get_p() {return p;}
  virtual void set_p(double _p){
    this->p=_p;
    geo=geoD(p);
  }
  virtual double beta() {
    double x1, x2;
    gam=gamD(a,1.); 
    x1=(gam)(gen);
    gam=gamD(b,1.); 
    x2=(gam)(gen);
    return x1/(x1+x2);
  }
  virtual void set_beta(double _a, double _b){this->a=_a; this->b=_b;}
  virtual double get_a(){return a;}
  virtual double get_b(){return b;}
  virtual std::vector<double> dirichlet(){
    size_t k=alpha.size();
    std::vector<double> draw(k);
    double gsum=0;
    for(size_t j=0;j<k;j++){
      gam=gamD(alpha[j],1.); 
      draw[j]=(gam)(gen);
      gsum+=draw[j];
    }
    for(size_t j=0;j<k;j++) draw[j]/=gsum;
    return draw;
  }
  virtual void set_alpha(std::vector<double>& _alpha){this->alpha=_alpha;}
  virtual std::vector<double> get_alpha(){return alpha;}
 private:
  double df,shape, p, a, b;
  unsigned int n1, n2;
  std::vector<double> wts, alpha;
  genD gen;
  norD nor;
  uniD uni;
  chiD chi;
  gamD gam;
//  disD dis;
  geoD geo;
};

#elif defined (RNG_Rmath) 

#include <Rmath.h>

//abstract random number generator based on Rmath
class arn: public rn
{
 public:
  //constructor
 arn():df(1) {}
 arn(unsigned int n1, unsigned int n2):df(1) {::set_seed(n1, n2);}
  //virtual
  virtual ~arn() {}
  virtual double normal() {return ::norm_rand();}
  virtual double uniform() { return ::unif_rand();}
  virtual double chi_square() {return ::rchisq(df);}
  virtual double exp() {return ::exp_rand();}
  virtual double gamma() {return ::rgamma(shape,1.);}
  virtual size_t discrete() {
    size_t p=wts.size(), x=0;
    std::vector<int> vOut (p,0);
    ::rmultinom(1,&wts[0],p,&vOut[0]); 
    if(vOut[0]==0) for(size_t j=1;j<p;j++) x += j*vOut[j]; 
    return x;
  }
/*
    size_t p=wts.size();
    std::vector<int> vOut (p,0);
    std::vector<int> seq (p,0);
    for(size_t j=0;j<p;j++) seq[j]=j; // creates vector of sequence 1,...,p
    ::rmultinom(1,&wts[0],p,&vOut[0]); // vector vOut contains selected component
    return std::inner_product(vOut.begin(),vOut.end(),seq.begin(),0);
  }
*/
  virtual size_t geometric() {return ::rgeom(p);}
  virtual void set_df(double _df) {this->df=_df;}
  virtual double get_df() {return df;}
  virtual void set_shape(double _shape) {this->shape=_shape;}
  virtual double get_shape() {return shape;}
  virtual void set_wts(std::vector<double>& _wts) {
    double smw=0.;
    wts.clear();
    for(size_t j=0;j<_wts.size();j++) smw+=_wts[j];
    for(size_t j=0;j<_wts.size();j++) wts.push_back(_wts[j]/smw);
  }
  virtual double get_p() {return p;}
  virtual void set_p(double _p){this->p=_p;}
  void set_seed(unsigned int n1, unsigned int n2) 
  {::set_seed(n1, n2);}
  void get_seed(unsigned int* n1, unsigned int* n2) 
  {::get_seed(n1, n2);}
  virtual double beta() {return ::rbeta(a, b);}
  virtual void set_beta(double _a, double _b){this->a=_a; this->b=_b;}
  virtual double get_a(){return a;}
  virtual double get_b(){return b;}
  virtual std::vector<double> dirichlet(){
    size_t k=alpha.size();
    std::vector<double> draw(k);
    double gsum=0;
    for(size_t j=0;j<k;j++){
      draw[j]=::rgamma(alpha[j],1.);
      gsum+=draw[j];
    }
    for(size_t j=0;j<k;j++) draw[j]/=gsum;
    return draw;
  }
  virtual void set_alpha(std::vector<double>& _alpha){this->alpha=_alpha;}
  virtual std::vector<double> get_alpha(){return alpha;}
 private:
  double df,shape,p,a,b;
  std::vector<double> wts,alpha;
};

#else // YesRcpp

//abstract random number generator based on R/Rcpp
class arn: public rn
{
 public:
  //constructor
 arn():df(1) {}
  //virtual
  virtual ~arn() {}
  virtual double normal() {return R::norm_rand();}
  virtual double uniform() { return R::unif_rand();}
  virtual double chi_square() {return R::rchisq(df);}
  virtual double exp() {return R::exp_rand();}
  virtual double gamma() {return R::rgamma(shape,1.);}
  virtual size_t discrete() {
    size_t p=wts.size(), x=0;
    std::vector<int> vOut (p,0);
    R::rmultinom(1,&wts[0],p,&vOut[0]); 
    if(vOut[0]==0) for(size_t j=1;j<p;j++) x += j*vOut[j]; 
    return x;
  }
/*
    size_t p=wts.size();
    std::vector<int> vOut (p,0); // vector of multionomial output
    std::vector<int> seq (p,0); 
    for(size_t j=0;j<p;j++) seq[j]=j;
    R::rmultinom(1,&wts[0],p,&vOut[0]); // vector vOut contains selected component
    return std::inner_product(vOut.begin(),vOut.end(),seq.begin(),0);
  }
*/
  virtual size_t geometric() {return R::rgeom(p);}
  virtual void set_df(double _df) {this->df=_df;}
  virtual double get_df() {return df;}
  virtual void set_shape(double _shape) {this->shape=_shape;}
  virtual double get_shape() {return shape;}
  virtual void set_wts(std::vector<double>& _wts) {
    double smw=0.;
    wts.clear();
    for(size_t j=0;j<_wts.size();j++) smw+=_wts[j];
    for(size_t j=0;j<_wts.size();j++) wts.push_back(_wts[j]/smw);
  }
  virtual double get_p() {return p;}
  virtual void set_p(double _p){this->p=_p;}
  virtual double beta() {return R::rbeta(a, b);}
  virtual void set_beta(double _a, double _b){this->a=_a; this->b=_b;}
  virtual double get_a(){return a;}
  virtual double get_b(){return b;}
  virtual std::vector<double> dirichlet(){
    size_t k=alpha.size();
    std::vector<double> draw(k);
    double gsum=0;
    for(size_t j=0;j<k;j++){
      draw[j]=R::rgamma(alpha[j],1.);
      gsum+=draw[j];
    }
    for(size_t j=0;j<k;j++) draw[j]/=gsum;
    return draw;
  }
  virtual void set_alpha(std::vector<double>& _alpha){this->alpha=_alpha;}
  virtual std::vector<double> get_alpha(){return alpha;}
 private:
  double df,shape,p,a,b;
  std::vector<double> wts,alpha;
  Rcpp::RNGScope RNGstate;
};

#endif 

#endif 
