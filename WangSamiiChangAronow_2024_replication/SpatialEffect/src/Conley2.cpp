// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export()]]
Rcpp::List ConleySE2(arma::colvec res, arma::mat W_meat, arma::mat dist_kernel, arma::mat XXinv, int kernel, int trim,int if_edof){

  // Define variables
  int N = res.size();
  // double mu=1;
  // double v=2;
  double muprime=1;
  double vprime=2;
  // arma::mat B = arma::zeros(N, N);
  arma::mat Bprime = arma::zeros(N, N);

  // arma::mat bread = inv(W_meat.t() * W_meat);
  arma::mat meat = arma::zeros(W_meat.n_cols, W_meat.n_cols);
  // arma::mat meat2 = arma::zeros(W_meat.n_cols, W_meat.n_cols);
  arma::mat meat_add = arma::zeros(W_meat.n_cols, W_meat.n_cols);
  //arma::mat meat_EHW = arma::zeros(W_meat.n_cols, W_meat.n_cols);
  arma::vec w = arma::zeros(2);
  arma::vec z_vec = arma::zeros(N);
  arma::mat Identity_matrix = arma::eye(N,N);
  arma::mat M_mat = arma::zeros(N,N);

    //coef vector
  w(0)=0;
  w(1)=1;

  // clock.tick("rmm");
  //Residual Maker Matrix
  M_mat = Identity_matrix - W_meat * XXinv * W_meat.t();

  z_vec = W_meat * XXinv * w;
 
  // clock.tock("rmm");
  
  // clock.tick("edof");
  // clock.tick("edof1_for");

if (if_edof==1){


for (int i=0;i<N;i++){
 // for (int j=0;j<N;j++){
  //    Bprime.col(i) += M_mat.col(j) * (z_vec(j) * dist_kernel(i,j) )/XXinv(1,1)*z_vec(i);

 // }



// cout << "mu: " << muprime << endl;
// cout << "v: " << vprime << endl;
//   clock.tick("percent");
//   Bprime.col(i) = M_mat * (z_vec % dist_kernel.col(i) )/XXinv(1,1)*z_vec(i);
//   clock.tock("percent");
// muprime= trace(Bprime);
// vprime= 2*trace(Bprime*Bprime);

// cout << "mu: " << muprime << endl;
// cout << "v: " << vprime << endl;
  // clock.tick("2n");
  Bprime.col(i) = (z_vec % dist_kernel.col(i) )/XXinv(1,1)*z_vec(i)-W_meat * XXinv * (W_meat.t() * (z_vec % dist_kernel.col(i) ))/XXinv(1,1)*z_vec(i);
  // clock.tock("2n");
}

  // clock.tock("edof1_for");
  
  // clock.tick("edof2_for");
muprime= trace(Bprime);
vprime= 2*trace(Bprime*Bprime);

// cout << "mu: " << muprime << endl;
// cout << "v: " << vprime << endl;
  // clock.tock("edof2_for");
}  

//   // Effective Degrees of Freedom
//   if (if_edof==1){
//     clock.tick("edof1");
//     B = M_mat * (dist_kernel % (z_vec * z_vec.t())) * M_mat/XXinv(1,1);
//     clock.tock("edof1");
//     clock.tick("edof2");
//     mu = trace(B);
//     clock.tock("edof2");
//     clock.tick("edof3");
//     v = 2*trace(B*B);
//     clock.tock("edof3");
//   }
//  clock.tock("edof");

//   clock.tick("conley1");

// cout << "mu: " << mu << endl;

// cout << "v: " << v << endl;

// cout << "muprime: " << muprime << endl;

// cout << "vprime: " << vprime << endl;  
  for (int i = 0; i < N; i++){
    meat_add = (res(i) * W_meat.row(i)).t() * (((res % (dist_kernel.col(i))).t()) * W_meat);
    meat = meat + meat_add;
    //meat_EHW = meat_EHW + pow(res(i), 2) * W_meat.row(i).t() * W_meat.row(i);
  }

  // clock.tock("conley1");


  // clock.tick("conley2");

  // meat2 = W_meat.t() * ((res* res.t() ) % dist_kernel) * W_meat;
  // // arma::mat VCE = bread.t() * meat * bread;
  // arma::mat VCE_EHW = bread.t() * meat_EHW * bread;

  // clock.tock("conley2");

  // // storage

    return Rcpp::List::create(Rcpp::Named("VCE_meat2") = meat, 
                              Rcpp::Named("mu") = muprime,
                              Rcpp::Named("v") = vprime);
 
}
