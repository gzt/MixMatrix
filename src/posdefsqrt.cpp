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
  if ((min(eigvalU) < 1e-7)) Rcpp::stop("error: possibly singular input");
  arma::eig_sym(eigvalV, eigvecV, V);
  if ((min(eigvalV) < 1e-7)) Rcpp::stop("error: possibly singular input");
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
  if ((min(eigvalU) < 1e-7)) Rcpp::stop("error: possibly singular input");
  arma::eig_sym(eigvalV, eigvecV, V);
  if ((min(eigvalV) < 1e-7)) Rcpp::stop("error: possibly singular input");
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
  arma::mat Uinv;
  bool res = arma::inv_sympd(Uinv, U);
  if(!res) stop("error: singular or non-positive definite input");
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
  arma::mat Uinv;
  bool res = arma::inv_sympd(Uinv, U);
  if(!res) stop("error: singular or non-positive definite input");
  arma::cube results(p, p, numslices);
  for(int i = 0; i < numslices; i++){
    results.slice(i) = trans(x.slice(i)) * Uinv * x.slice(i);
  }
  return results;
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
arma::colvec dmat_inv_t_calc(arma::cube & x, double df, arma::mat & mean,
                         arma::mat & U, arma::mat & V){
  // note this needs to be scaled by a constant outside of this function!
  arma::colvec eigvalU;
  arma::mat eigvecU;
  arma::colvec eigvalV;
  arma::mat eigvecV;
  arma::eig_sym(eigvalU, eigvecU, U);
  if ((min(eigvalU) < 1e-7)) Rcpp::stop("error: possibly singular input");
  arma::eig_sym(eigvalV, eigvecV, V);
  if ((min(eigvalV) < 1e-7)) Rcpp::stop("error: possibly singular input");
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
    arma::mat m = identitymat - Uinv * XM * Vinv * trans(XM) ;
    double mval, sign;
    log_det(mval, sign, m);
    if(sign <= 0){
      logresult[i] = R_NaN;
      warning("warning: probability distribution undefined when det < 0. observation: %d ", i+1);
    } else{
    logresult[i] = -.5 * p * logdetU - .5 * n * logdetV - .5 * (df - 2) * mval;
    }
  }
  return logresult;
}


// [[Rcpp::export]]
arma::cube cubeinv(arma::cube & x){
  // generates a cube of X_i^{-1} for i = 1:n, X_i is p x p
  int p = x.n_cols;
  int numslices = x.n_slices;
  if (x.n_rows != x.n_cols) Rcpp::stop("error: non-conformable dimensions");
  arma::cube results(p, p, numslices);
  for(int i = 0; i < numslices; i++){
    arma::mat Xinv;
    bool res = arma::inv_sympd(Xinv, x.slice(i));
    if(!res) stop("error: singular or non-positive definite input");
    results.slice(i) = Xinv;
  }
  return results;
}

// [[Rcpp::export]]
arma::cube cubemult(arma::cube & x, arma::cube & y){
  // multiplies t(x) * y by slice
    int q  = x.n_cols;
  int r = y.n_cols;
  int numslices = x.n_slices;
  if (y.n_rows != x.n_rows) Rcpp::stop("error: non-conformable dimensions");
  arma::cube results(q, r, numslices);
  for(int i = 0; i < numslices; i++){
    results.slice(i) = trans(x.slice(i)) * y.slice(i);
  }
  return results;
}

// [[Rcpp::export]]
double detsum(arma::cube & x){
  // takes log det of each slice and adds it up
  if (x.n_rows != x.n_cols) Rcpp::stop("error: non-conformable dimensions");
  int numslices = x.n_slices;
  double result = 0;
  double mval, sign;
  for(int i = 0; i < numslices; i++){
  log_det(mval, sign, x.slice(i));
  if(sign <= 0){
      Rcpp::stop("error: result undefined when det < 0. observation: %d ", i+1);
    }
  result += mval;
  }
  return result;
}


// [[Rcpp::export]]
arma::cube axbt(arma::cube & a, arma::mat & x, arma::cube & b){
  // a is p x q x n, x is q x q, b is r x q x n; returns a * x * t(b) cube
  int p = a.n_rows;
  // int q  = x.n_cols;
  int r = b.n_rows;
  int numslices = a.n_slices;
  if (x.n_rows != x.n_cols || a.n_cols != x.n_rows || x.n_cols != b.n_cols) Rcpp::stop("error: non-conformable dimensions");
  arma::cube results(p, r, numslices);
  for(int i = 0; i < numslices; i++){
    results.slice(i) = (a.slice(i)) * x * trans(b.slice(i));
  }
  return results;
}

