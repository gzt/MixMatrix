#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat posdefsqrt(arma::mat & x){
  arma::mat eigmat;
  arma::vec eigvec;
  arma::eig_sym(eigvec, eigmat, x);
  if(min(eigvec) < 1e-7) throw Rcpp::exception("error: possibly singular matrix");
  return eigmat * diagmat(sqrt(eigvec)) * eigmat.t();
}


// [[Rcpp::export]]
arma::mat posdefinvsqrt(arma::mat & x){
  arma::mat eigmat;
  arma::vec eigvec;
  arma::eig_sym(eigvec, eigmat, x);
  if(min(eigvec) < 1e-7) throw Rcpp::exception("error: possibly singular matrix");
  return eigmat * diagmat(1/sqrt(eigvec)) * eigmat.t();
}

// [[Rcpp::export]]
arma::cube rmat_inv_t_calc(arma::cube & S, arma::cube & mat,
                            arma::mat & U, arma::mat & V,
                            arma::mat & mean){
  //S <- stats::rWishart(n, df + dims[1] - 1, diag(dims[1]))
  int p = S.n_rows;
  int q = V.n_rows;
  int n = S.n_slices;

  arma::mat Usqrt = posdefsqrt(U);
  arma::mat Vsqrt = posdefsqrt(V);

  arma::cube SXX(p,p,n);

  for (int i = 0; i < n; i++) {
    SXX.slice(i) = (S.slice(i) + mat.slice(i) * trans(mat.slice(i)));
    SXX.slice(i) = posdefinvsqrt(SXX.slice(i));
  }
  arma::cube result(p,q,n);
    for (int i = 0; i < n; i++) {
      result.slice(i) = Usqrt * SXX.slice(i)* mat.slice(i) * Vsqrt + mean;
    }
    return result;
}



// [[Rcpp::export]]
arma::cube cubeinv(arma::cube & x){
  // generates a cube of X_i^{-1} for i = 1:n, X_i is p x p
  int p = x.n_cols;
  int numslices = x.n_slices;
  if (x.n_rows != x.n_cols) throw Rcpp::exception("error: non-conformable dimensions");
  arma::cube results(p, p, numslices);
  for(int i = 0; i < numslices; i++){
    arma::mat Xinv;
    bool res = arma::inv_sympd(Xinv, x.slice(i));
    if(!res) throw Rcpp::exception("error: singular or non-positive definite input");
    results.slice(i) = Xinv;
  }
  return results;
}


// [[Rcpp::export]]
arma::cube xatx(arma::cube & x, arma::mat & U){
  // generates X * solve(U) * t(X)
  int n = x.n_rows;
  int numslices = x.n_slices;
  if (x.n_cols != U.n_rows || x.n_cols != U.n_cols) throw Rcpp::exception("error: non-conformable dimensions");
  arma::mat Uinv;
  bool res = arma::inv_sympd(Uinv, U);
  if(!res) throw Rcpp::exception("error: singular or non-positive definite input");
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
  if (x.n_rows != U.n_rows || x.n_rows != U.n_cols) throw Rcpp::exception("error: non-conformable dimensions");
  arma::mat Uinv;
  bool res = arma::inv_sympd(Uinv, U);
  if(!res) throw Rcpp::exception("error: singular or non-positive definite input");
  arma::cube results(p, p, numslices);
  for(int i = 0; i < numslices; i++){
    results.slice(i) = trans(x.slice(i)) * Uinv * x.slice(i);
  }
  return results;
}


