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
bool testsymmetric(arma::mat & x, double tol){
  int nrow = x.n_rows;
  int ncol = x.n_cols;
  double total = 0.0;
  arma::mat flop = abs((x - x.t()));
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
     total += (flop(i,j));
    }
  }
  return (total < tol);
}

// [[Rcpp::export]]
arma::mat txax(arma::mat & x, arma::mat & A){
  // outputs t(x) %*% A %*% x
  int nrow = x.n_rows;
  // int ncol = x.n_cols;
  int nArow = A.n_rows;
  int nAcol = A.n_cols;
  if(nArow != nAcol || nrow != nArow) stop("error: non-conformable dimensions");
  return(x.t() * A * x);
}

// [[Rcpp::export]]
arma::mat xatx(arma::mat & x, arma::mat & A){
  // outputs (x) %*% A %*% t(x)
  // int nrow = x.n_rows;
  int ncol = x.n_cols;
  int nArow = A.n_rows;
  int nAcol = A.n_cols;
  if(nArow != nAcol || ncol != nArow) stop("error: non-conformable dimensions");
  return(x * A * x.t());
}

