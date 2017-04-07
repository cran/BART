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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using std::endl;

#ifdef _OPENMP
#include <omp.h>
#endif

//pure virtual base class for random numbers
class rn
{
public:
   virtual double normal() = 0; //standard normal
   virtual double uniform() = 0; //uniform(0,1)
   virtual double chi_square() = 0; //chi-square
   virtual double exp() = 0; //exponential
   virtual void set_df(int df) = 0; //set df for chi-square
   virtual ~rn() {}
};

#ifdef NoRcpp

#include <random> 

using std::cout;

//abstract random number generator based on C++ <random>
class arn: public rn
{
  //typedefs
  typedef std::default_random_engine genD;
  typedef std::normal_distribution<double> norD;
  typedef std::uniform_real_distribution<double> uniD;
  typedef std::chi_squared_distribution<double> chiD;
 public:
  //constructor
  arn();
  //virtual
  virtual ~arn();
  virtual double normal() {return (*nor)(*gen);}
  virtual double uniform() {return (*uni)(*gen);}
  virtual double chi_square() {return (*chi)(*gen);}
  virtual double exp() {return -log(uniform());}
  virtual void set_df(int df);
  int get_df()  {return df;}
  void set_seed(int seed);
 private:
  int df;
  genD* gen;
  norD* nor;
  uniD* uni;
  chiD* chi;
};

#define PI 3.1415926535897931

#else // YesRcpp

#include <Rcpp.h>

#define printf Rprintf
#define cout Rcpp::Rcout

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
  virtual double chi_square() {return R::rchisq((double)df);}
  virtual double exp() {return R::exp_rand();}
  virtual void set_df(int df) {this->df=df;}
  int get_df() {return df;}
 private:
  int df;
  Rcpp::RNGScope RNGstate;
};

#endif 

// log(2*pi)
#define LTPI 1.83787706640934536

#endif 
