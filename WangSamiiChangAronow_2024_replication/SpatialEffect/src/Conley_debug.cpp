#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export()]]
void ConleySE(arma::colvec res, arma::mat W_meat, arma::mat dist, arma::mat XXinv, double cutoff, int kernel, int trim,int if_edof){

  // Define variables
  int N = res.size();
  double u = 0;
  double mu=1;
  double v=2;
  arma::mat B = arma::zeros(N, N);
  // arma::mat bread = inv(W_meat.t() * W_meat);
  arma::mat meat = arma::zeros(W_meat.n_cols, W_meat.n_cols);
  arma::mat meat_add = arma::zeros(W_meat.n_cols, W_meat.n_cols);
  arma::mat meat_EHW = arma::zeros(W_meat.n_cols, W_meat.n_cols);
  arma::mat dist_kernel = arma::zeros(N, N);
  arma::vec w = arma::zeros(2);
  arma::vec z_vec = arma::zeros(N);
  arma::mat Identity_matrix = arma::eye(N,N);
  arma::mat M_mat = arma::zeros(N,N);

    //coef vector
  w(0)=0;
  w(1)=1;
  //Residual Maker Matrix
  M_mat = Identity_matrix - W_meat * XXinv * W_meat.t();
  z_vec = W_meat * XXinv * w;

  //Calculate Kernel Matrix
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      u = dist(i, j)/cutoff;

      
      if (kernel == 3) {
        dist_kernel(i, j) = (0.75*(1-pow(u, 2))*(u <= 1));
      }
      else if (kernel == 1) {
        dist_kernel(i, j) = (u <= 1);
      }
      else {
        dist_kernel(i, j) = ((1-u)*(u <= 1));
      }
    }

  }


  
}
