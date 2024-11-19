#include <RcppArmadillo.h>
#include <random>


// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;
// [[Rcpp::export()]]

NumericVector nonmonoeffect_c(NumericVector d, double sh, double sc, double a){
  
  
  return ( (Rcpp::dgamma(d, sh, 2*sc,0) - a*Rcpp::dgamma(d, 5*sh, sc,0))*((-1)*d*d/36+1));
}
// [[Rcpp::export]]


Rcpp::List effect_nonmo1(arma::colvec Y0, NumericMatrix dist, arma::mat Zdata, arma::colvec alpha, double effect,double sh, double sc, double a){

    // Define variables
    int Ny = Y0.size();
    int Nsim = Zdata.n_cols;
    int Nz = Zdata.n_rows;

    // Final outoput matrix
    arma::mat output(Ny, Nsim);
    NumericMatrix effect_by_distance(Ny,Nz);

    for (int i=0; i<Ny; i++){
            effect_by_distance(i,_) = nonmonoeffect_c(dist(i,_), sh, sc, a);
        }
    
    arma::mat effect_by_distance_mat = Rcpp::as<arma::mat>(effect_by_distance);
    
    output = effect_by_distance_mat * Zdata;
    for ( int i=0; i<Ny ; i++){

        for (int j=0; j<Nsim; j++){

            output(i,j) =Y0(i) + alpha(i) * effect*output(i,j);
           
        }
    }
    
    // storage
    return Rcpp::List::create(Rcpp::Named("Sim_data") = output);



}


