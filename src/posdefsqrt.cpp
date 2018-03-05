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



// [[Rcpp::export]]
arma::colvec dmatnorm_calc(arma::cube & x, arma::mat & mean,
                           arma::mat & U, arma::mat & V){
  arma::colvec eigvalU;
  arma::mat eigvecU;
  arma::colvec eigvalV;
  arma::mat eigvecV;
  arma::eig_sym(eigvalU, eigvecU, U);
  if ((min(eigvalU) < 1e-7)) Rcpp::stop("error: possibly non-singular input");
  arma::eig_sym(eigvalV, eigvecV, V);
  if ((min(eigvalV) < 1e-7)) Rcpp::stop("error: possibly non-singular input");
  int n = x.n_rows;
  int p = x.n_cols;
  int numslices = x.n_slices;
  // char* uchar = "U";
  double logdetU = 2*sum(log(eigvalU));
  double logdetV = 2*sum(log(eigvalV));
  arma::mat Uinv = eigvecU * diagmat(1/eigvalU) * eigvecU.t();
  arma::mat Vinv =  eigvecV * diagmat(1/eigvalV) * eigvecV.t();
  arma::cube XM(x.n_rows, x.n_cols, x.n_slices);
  for (int i = 0; i < numslices; i++) {
    XM.slice(i) = x.slice(i) - mean;
  }
  arma::colvec logresult(x.n_slices);
  for (int i = 0; i < numslices; i++) {
    logresult(i) = -0.5 * n * p * log(2 * PI) - 0.5 * n * logdetU -
      0.5 * p * logdetV - 0.5 * trace( Vinv * trans(XM.slice(i)) * Uinv * XM.slice(i));
  }
  return logresult;
}

