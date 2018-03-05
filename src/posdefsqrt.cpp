#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat posdefsqrt(arma::mat & x){
  arma::mat eigmat;
  arma::vec eigvec;

  arma::eig_sym(eigvec, eigmat, x);
  if(min(eigvec) < 1e-7) stop("error: possibly singular matrix");
  return eigmat * diagmat(sqrt(eigvec)) * eigmat.t();
}

// [[Rcpp::export]]
arma::mat posdefinvsqrt(arma::mat & x){
  arma::mat eigmat;
  arma::vec eigvec;

  arma::eig_sym(eigvec, eigmat, x);
  if(min(eigvec) < 1e-7) stop("error: possibly singular matrix");
  return eigmat * diagmat(1/sqrt(eigvec)) * eigmat.t();
}

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
