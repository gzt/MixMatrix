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
  double logdetU = sum(log(eigvalU));
  double logdetV = sum(log(eigvalV));
  arma::mat Uinv = eigvecU * diagmat(1/eigvalU) * eigvecU.t();
  arma::mat Vinv = eigvecV * diagmat(1/eigvalV) * eigvecV.t();
  arma::colvec logresult(x.n_slices);
  for (int i = 0; i < numslices; i++) {
    arma::mat XM = x.slice(i) - mean;
    logresult(i) = -0.5 * n * p * log(2 * PI) - 0.5 * n * logdetU -
      0.5 * p * logdetV - 0.5 * trace( Vinv * trans(XM) * Uinv * XM);
  }
  return logresult;
}

// [[Rcpp::export]]
arma::colvec dmat_t_calc(arma::cube & x, double df, arma::mat & mean,
                           arma::mat & U, arma::mat & V){
  // note this needs to be scaled by a constant outside of this function!
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

  double logdetU = sum(log(eigvalU));
  double logdetV = sum(log(eigvalV));
  arma::mat Uinv = eigvecU * diagmat(1/eigvalU) * eigvecU.t();
  arma::mat Vinv = eigvecV * diagmat(1/eigvalV) * eigvecV.t();
  arma::colvec logresult(x.n_slices);

  for (int i = 0; i < numslices; i++) {
    arma::mat XM = x.slice(i) - mean;
    arma::mat identitymat = arma::eye(n, n);
    arma::mat m = identitymat + Uinv * XM * Vinv * trans(XM) ;
    logresult[i] = -.5 * p * logdetU - .5 * n * logdetV - .5 * (df + n + p -1) * log(det(m));
  }

  return logresult;

}

// [[Rcpp::export]]
arma::cube xatx(arma::cube & x, arma::mat & U){
  // generates X * solve(U) * t(X)
  int n = x.n_rows;
  int numslices = x.n_slices;
  if (x.n_cols != U.n_rows || x.n_cols != U.n_cols) Rcpp::stop("error: non-conformable dimensions");
  arma::mat Uinv = arma::inv_sympd(U);
  arma::cube results(n,n,numslices);
  for(int i = 0; i < numslices; i++){
    results.slice(i) = (x.slice(i)) * Uinv * trans(x.slice(i));
  }
  return results;
}


// [[Rcpp::export]]
arma::cube txax(arma::cube & x, arma::mat & U){
  // generates t(X) * solve(U) * X
  int p = x.n_cols;
  int numslices = x.n_slices;
  if (x.n_rows != U.n_rows || x.n_rows != U.n_cols) Rcpp::stop("error: non-conformable dimensions");
  arma::mat Uinv = arma::inv_sympd(U);
  arma::cube results(p,p,numslices);
  for(int i = 0; i < numslices; i++){
    results.slice(i) = trans(x.slice(i)) * Uinv * x.slice(i);
  }
  return results;
}
