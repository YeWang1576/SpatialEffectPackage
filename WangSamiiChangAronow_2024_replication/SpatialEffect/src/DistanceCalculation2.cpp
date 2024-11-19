#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export()]]
Rcpp::List DistanceCalculation2(arma::colvec x_coord, arma::colvec y_coord, arma::colvec x_coord2, arma::colvec y_coord2, int metric){

  // Define variables
  int Ny = x_coord.size();
  int Nz = x_coord2.size();
  
  arma::mat outcome_coords = join_rows(x_coord, y_coord);
  arma::mat node_coords = join_rows(x_coord2, y_coord2);
  
  arma::mat dist = arma::zeros(Ny, Nz);  

  for (int i = 0; i < Ny; i++){
    arma::colvec point1 = outcome_coords.row(i).t();
    for (int j = 0; j < Nz; j++){
      arma::colvec point2 = node_coords.row(j).t();
      if (metric == 1){
        dist(i, j) = sqrt(sum(pow(point1 - point2, 2)));
      }
    }
  }

  // storage
  return Rcpp::List::create(Rcpp::Named("Dist_mat") = dist);
}
