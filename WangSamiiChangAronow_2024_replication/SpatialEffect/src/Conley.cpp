#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export()]]
Rcpp::List ConleySE(arma::colvec res, arma::mat W_meat, arma::mat dist, arma::mat XXinv, double cutoff, int kernel, int trim,int if_edof){

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


  //Trimming
   if (trim == 1){

      arma::mat U;
      arma::vec s;
      arma::eig_sym( s, U, dist_kernel );
//      cout << s << endl;
     // cout << U<<endl;
      for (int k = 0; k < s.n_elem; k++){
      if (s(k) < 0){
         s(k) = 0;
      }
     } 
     arma::mat S = diagmat(s);
     dist_kernel = U * S * U.t();
   }

  // Effective Degrees of Freedom
  if (if_edof==1){

    B = M_mat * (dist_kernel % (z_vec * z_vec.t())) * M_mat/XXinv(1,1);
    mu = trace(B);
    v = 2*trace(B*B);
  }

  for (int i = 0; i < N; i++){
    meat_add = (res(i) * W_meat.row(i)).t() * (((res % (dist_kernel.col(i))).t()) * W_meat);
    meat = meat + meat_add;
    meat_EHW = meat_EHW + pow(res(i), 2) * W_meat.row(i).t() * W_meat.row(i);
  }
  
  // arma::mat VCE = bread.t() * meat * bread;
  // arma::mat VCE_EHW = bread.t() * meat_EHW * bread;

  // storage
  if (cutoff == 0){
    return Rcpp::List::create(Rcpp::Named("VCE_meat") = meat_EHW, 
                              Rcpp::Named("Dist_kernel") = dist_kernel,
                              Rcpp::Named("mu") = mu,
                              Rcpp::Named("v") = v);
  }else{
    return Rcpp::List::create(Rcpp::Named("VCE_meat") = meat, 
                              Rcpp::Named("Dist_kernel") = dist_kernel,
                              Rcpp::Named("mu") = mu,
                              Rcpp::Named("v") = v);
  }
}
