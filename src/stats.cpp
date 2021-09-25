
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
arma::colvec dmat_inv_t_calc(arma::cube& x, double df, arma::mat& mean,
                             arma::mat& U, arma::mat& V) {
  // note this needs to be scaled by a constant outside of this function!
  arma::colvec eigvalU;
  arma::mat eigvecU;
  arma::colvec eigvalV;
  arma::mat eigvecV;
  arma::eig_sym(eigvalU, eigvecU, U);
  if ((min(eigvalU) < 1e-7))
    throw Rcpp::exception("error: possibly singular input");
  arma::eig_sym(eigvalV, eigvecV, V);
  if ((min(eigvalV) < 1e-7))
    throw Rcpp::exception("error: possibly singular input");
  int n = x.n_rows;
  int p = x.n_cols;
  int numslices = x.n_slices;
  double logdetU = sum(log(eigvalU));
  double logdetV = sum(log(eigvalV));
  arma::mat Uinv = eigvecU * diagmat(1 / eigvalU) * eigvecU.t();
  arma::mat Vinv = eigvecV * diagmat(1 / eigvalV) * eigvecV.t();
  arma::colvec logresult(x.n_slices);
  for (int i = 0; i < numslices; i++) {
    arma::mat XM = x.slice(i) - mean;
    arma::mat identitymat = arma::eye(n, n);
    arma::mat m = identitymat - Uinv * XM * Vinv * trans(XM);
    double mval, sign;
    log_det(mval, sign, m);
    if (sign <= 0) {
      logresult[i] = R_NaN;
      // Rcpp::warning("warning: probability distribution undefined when det <
      // 0. observation: %d ", i+1);
    } else {
      logresult[i] =
          -.5 * n * logdetU - .5 * p * logdetV - .5 * (df - 2) * mval;
    }
  }
  return logresult;
}

// [[Rcpp::export]]
arma::colvec dmatnorm_calc(arma::cube& x, arma::mat& mean, arma::mat& U,
                           arma::mat& V) {
  arma::colvec eigvalU;
  arma::mat eigvecU;
  arma::colvec eigvalV;
  arma::mat eigvecV;
  arma::eig_sym(eigvalU, eigvecU, U);
  if ((min(eigvalU) < 1e-7))
    throw Rcpp::exception("error: possibly singular input");
  arma::eig_sym(eigvalV, eigvecV, V);
  if ((min(eigvalV) < 1e-7))
    throw Rcpp::exception("error: possibly singular input");
  int n = x.n_rows;
  int p = x.n_cols;
  int numslices = x.n_slices;
  // char* uchar = "U";
  double logdetU = sum(log(eigvalU));
  double logdetV = sum(log(eigvalV));
  arma::mat Uinv = eigvecU * diagmat(1 / eigvalU) * eigvecU.t();
  arma::mat Vinv = eigvecV * diagmat(1 / eigvalV) * eigvecV.t();
  arma::colvec logresult(x.n_slices);
  for (int i = 0; i < numslices; i++) {
    arma::mat XM = x.slice(i) - mean;
    logresult(i) = -0.5 * n * p * log(2 * M_PI) - 0.5 * p * logdetU -
                   0.5 * n * logdetV -
                   0.5 * trace(Vinv * trans(XM) * Uinv * XM);
  }
  return logresult;
}

// [[Rcpp::export]]
arma::colvec dmat_t_calc(arma::cube& x, double df, arma::mat& mean,
                         arma::mat& U, arma::mat& V) {
  // note this needs to be scaled by a constant outside of this function!
  arma::colvec eigvalU;
  arma::mat eigvecU;
  arma::colvec eigvalV;
  arma::mat eigvecV;
  arma::eig_sym(eigvalU, eigvecU, U);
  if ((min(eigvalU) < 1e-7))
    throw Rcpp::exception("error: possibly singular input");
  arma::eig_sym(eigvalV, eigvecV, V);
  if ((min(eigvalV) < 1e-7))
    throw Rcpp::exception("error: possibly singular input");
  int n = x.n_rows;
  int p = x.n_cols;
  int numslices = x.n_slices;
  double logdetU = sum(log(eigvalU));
  double logdetV = sum(log(eigvalV));
  arma::mat Uinv = eigvecU * diagmat(1 / eigvalU) * eigvecU.t();
  arma::mat Vinv = eigvecV * diagmat(1 / eigvalV) * eigvecV.t();
  arma::colvec logresult(x.n_slices);
  for (int i = 0; i < numslices; i++) {
    arma::mat XM = x.slice(i) - mean;
    arma::mat identitymat = arma::eye(n, n);
    arma::mat m = identitymat + Uinv * XM * Vinv * trans(XM);
    logresult[i] = -.5 * p * logdetU - .5 * n * logdetV -
                   .5 * (df + n + p - 1) * log(det(m));
  }
  return logresult;
}
