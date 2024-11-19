#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export()]]
Rcpp::List DistanceCalculation(arma::colvec x_coord, arma::colvec y_coord, int metric){

  // Define variables
  int N = x_coord.size();
  arma::mat coords = join_rows(x_coord, y_coord);
  arma::mat dist = arma::zeros(N, N);  

  for (int i = 0; i < N; i++){
    arma::colvec point1 = coords.row(i).t();
    for (int j = 0; j < N; j++){
      arma::colvec point2 = coords.row(j).t();
      if (metric == 1){
        dist(i, j) = sqrt(sum(pow(point1 - point2, 2)));
      }else if (metric == 2){
        Environment ddd("package:geosphere"); 
        Function distGeoFunc = ddd["distGeo"];
        NumericVector geo_distance = distGeoFunc(_["p1"] = point1.t(), _["p2"] = point2.t());
        dist(i, j) = geo_distance[0];
      }
    }
  }

  // storage
  return Rcpp::List::create(Rcpp::Named("Dist_mat") = dist);
}
