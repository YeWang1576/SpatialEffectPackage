#include <RcppArmadillo.h>
#include <random>


// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;
// [[Rcpp::export()]]

NumericVector nonmonoeffect_c(NumericVector d, double sh, double sc, double a){

    
    NumericVector weights(d.size(),1);

    for (int i = 0; i < d.size();i++) {
          
          if (d(i) > 6) {
            weights(i) = 0;
          }
    }

    NumericVector output = (Rcpp::dgamma(d, sh, 2*sc,0) - a*Rcpp::dgamma(d, 5*sh, sc,0))*((-1)*d*d/36+1);

    return ( output*weights);

    
}


//If nearest neighbor is treated, add an additinoal effect function
NumericVector effect2(NumericVector d, double sh, double sc, double a){

    NumericVector weights(d.size(),1);

    for (int i = 0; i < d.size(); i++) {
          
          if (d(i) > 6) {
            weights(i) = 0;
          }
    }
    NumericVector output = ( Rcpp::dgamma(d, 5*sh, sc,0))* ((-1)*d*d/36+1);
  //return ( Rcpp::dgamma(d, 5*sh, sc,0))*((-1)*d*d/36+1);
  return (output*weights);
   // NumericVector weights2(d.size(),0); 
   //return(weights2);
}


// [[Rcpp::export]]



Rcpp::List effect_interactive(arma::vec Y0, NumericMatrix dist, arma::mat Zdata, arma::ivec Zneighbor,arma::vec alpha, double effect,double sh, double sc, double a){

    // Define variables
    int Ny = Y0.size();
    int Nsim = Zdata.n_cols;
    int Nz = Zdata.n_rows;


    
    // Final outoput matrix
    arma::mat output(Ny, Nsim);
    arma::mat Zinteract(Nz, Nsim);



    //create interactive treatment effect indicator

    for (int j=0; j<Nsim; j++){

        for (int i=0; i<Nz; i++){
           
           Zinteract(i,j) = Zdata(i,j) * Zdata(Zneighbor(i)-1, j);
           //Zinteract(i,j) = Zdata(i,j) * Zdata(1, j);
        }
    }
    NumericMatrix effect_by_distance_c(Ny,Nz);
    NumericMatrix effect_by_distance_t(Ny,Nz);

    
     
    // If nearest neighbor is not treated
    for (int i=0; i<Ny; i++){
            effect_by_distance_c(i,_) = nonmonoeffect_c(dist(i,_), sh, sc, 1);
        }
    
    // If nearest neighbor is treated
    for (int i=0; i<Ny; i++){
            effect_by_distance_t(i,_) = effect2(dist(i,_), sh, sc, 0);
        }

    arma::mat effect_by_distance_mat_c = Rcpp::as<arma::mat>(effect_by_distance_c);
    arma::mat effect_by_distance_mat_t = Rcpp::as<arma::mat>(effect_by_distance_t);

    output = effect_by_distance_mat_c * Zdata + effect_by_distance_mat_t * Zinteract;
    
    for ( int i=0; i<Ny ; i++){

        for (int j=0; j<Nsim; j++){

            output(i,j) =Y0(i) + alpha(i) * effect*output(i,j);
           
        }
    }
    
    // storage
    return Rcpp::List::create(Rcpp::Named("Sim_data") = output);



}


