#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

// [[Rcpp::export()]]
Rcpp::List ConleySE_kernel(arma::mat dist, int N, double cutoff, int kernel, int trim){

// arma::mat bread = inv(W_meat.t() * W_meat);
  arma::mat dist_kernel = arma::zeros(N, N);
  double u = 0;
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

  
  if (cutoff == 0){
    return Rcpp::List::create(Rcpp::Named("Dist_kernel") = dist_kernel);
  }else{
    return Rcpp::List::create(Rcpp::Named("Dist_kernel") = dist_kernel);
  }
}
