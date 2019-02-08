#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// [[Rcpp::export()]]
arma::mat ConleySE(arma::colvec res, arma::mat X_mat, arma::colvec x_coord,
                arma::colvec y_coord, double cutoff, int metric){

  // Define variables
  int N = res.size();
  arma::mat bread = inv(X_mat.t() * X_mat);
  arma::mat meat = arma::zeros(X_mat.n_cols, X_mat.n_cols);
  arma::mat meat_EHW = arma::zeros(X_mat.n_cols, X_mat.n_cols);
  arma::mat coords = join_rows(x_coord, y_coord);

  for (int i = 0; i < N; i++){
    arma::colvec point1 = coords.row(i).t();
    arma::colvec dist = arma::zeros(N);
    for (int j = 0; j < N; j++){
      arma::colvec point2 = coords.row(j).t();
      if (metric == 1){
        dist(j) = sum(pow(point1 - point2, 2));
      }else if (metric == 2){
        Environment ddd("package:geosphere"); 
        Function distGeoFunc = ddd["distGeo"];
        NumericVector geo_distance = distGeoFunc(_["p1"] = point1.t(), _["p2"] = point2.t());
        dist(j) = geo_distance[0];
      }
    }
    //dist = sign(sign(sqrt(dist) - cutoff) + 1);
    arma::uvec dist_u = (sqrt(dist) < cutoff);
    meat = meat + (res(i) * X_mat.row(i)).t() * (((res % dist_u).t()) * X_mat);
    meat_EHW = meat_EHW + res(i) * X_mat.row(i).t() * X_mat.row(i);
  }
  //meat = meat / N;
  arma::mat VCE = bread.t() * meat * bread;
  arma::mat VCE_EHW = bread.t() * meat_EHW * bread;
  //arma::colvec Vdiag = VCE.diag();
  //arma::colvec Vdiag_EHW = VCE_EHW.diag();
  //for (int k = 0; k < Vdiag.size(); k++){
  //  if (Vdiag(k) <= 0){
  //    Vdiag(k) = Vdiag_EHW(k);
  //  }
  //}
  //arma::colvec ConleySE = sqrt(Vdiag);
  //arma::colvec EHWSE = sqrt(Vdiag_EHW);
  // storage
  if (cutoff == 0){
    return VCE_EHW;
  }else{
    return VCE;
  }
}
