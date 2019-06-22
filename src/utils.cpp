
/*  MixMatrix: Classification with Matrix Variate Normal and t distributions
 *  Copyright (C) 2018-9  GZ Thompson <gzthompson@gmail.com>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;



// [[Rcpp::export]]
bool testsymmetric(arma::mat x, double tol){
  int nrow = x.n_rows;
  int ncol = x.n_cols;
  double total = 0.0;
  arma::mat flop = abs((x - trans(x)));
  for(int i = 0; i < nrow; i++){
     for(int j = 0; j < ncol; j++){
      total += 1.0*(flop(i,j) > tol);
      if ( total >= 1) break;
     }
   }
  return (total < 1);
}


// [[Rcpp::export]]
arma::cube axbt(arma::cube & a, arma::mat & x, arma::cube & b){
  // a is p x q x n, x is q x q, b is r x q x n; returns a * x * t(b) cube
  int p = a.n_rows;
  // int q  = x.n_cols;
  int r = b.n_rows;
  int numslices = a.n_slices;
  if (x.n_rows != x.n_cols || a.n_cols != x.n_rows || x.n_cols != b.n_cols) throw Rcpp::exception("error: non-conformable dimensions");
  arma::cube results(p, r, numslices);
  for(int i = 0; i < numslices; i++){
    results.slice(i) = (a.slice(i)) * x * trans(b.slice(i));
  }
  return results;
}


// [[Rcpp::export]]
arma::cube cubemult(arma::cube & x, arma::cube & y){
  // multiplies t(x) * y by slice
    int q  = x.n_cols;
  int r = y.n_cols;
  int numslices = x.n_slices;
  if (y.n_rows != x.n_rows) throw Rcpp::exception("error: non-conformable dimensions");
  arma::cube results(q, r, numslices);
  for(int i = 0; i < numslices; i++){
    results.slice(i) = trans(x.slice(i)) * y.slice(i);
  }
  return results;
}

// [[Rcpp::export]]
double detsum(arma::cube & x){
  // takes log det of each slice and adds it up
  if (x.n_rows != x.n_cols) throw Rcpp::exception("error: non-conformable dimensions");
  int numslices = x.n_slices;
  double result = 0;
  double mval, sign;
  for(int i = 0; i < numslices; i++){
  log_det(mval, sign, x.slice(i));
  if(sign <= 0){
    throw Rcpp::exception("error: result undefined when det < 0. observation: %d ", i+1);
    }
  result += mval;
  }
  return result;
}
