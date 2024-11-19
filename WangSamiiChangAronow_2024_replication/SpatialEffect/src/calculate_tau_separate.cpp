#include <RcppArmadillo.h>
#include <random>


// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// [[Rcpp::export()]]

Rcpp::List calculate_tau_separate(arma::mat Y, arma::mat Zdata, arma::ivec Zneighbor){

    // Define variables
    int Ny = Y.n_rows;
    int Nz = Zdata.n_rows;

    // Final outoput matrix
    arma::mat output_t(Ny, Nz); //if neighbor is treated
    arma::mat output_c(Ny, Nz); //if neighbor is not treated

    
    for ( int i=0; i<Ny ; i++){

        for (int j=0; j<Nz; j++){

          
           output_t(i,j) = arma::sum(Y.row(i) % Zdata.row(j) % Zdata.row(Zneighbor(j)-1)) / arma::sum (Zdata.row(j) % Zdata.row(Zneighbor(j)-1)) - arma::sum(Y.row(i) % (1-Zdata.row(j))% Zdata.row(Zneighbor(j)-1))/ arma::sum ((1-Zdata.row(j))% Zdata.row(Zneighbor(j)-1));
           output_c(i,j) = arma::sum(Y.row(i) % Zdata.row(j) % (1-Zdata.row(Zneighbor(j)-1))) / arma::sum (Zdata.row(j) % (1-Zdata.row(Zneighbor(j)-1))) - arma::sum(Y.row(i) % (1-Zdata.row(j))% (1-Zdata.row(Zneighbor(j)-1)))/ arma::sum ((1-Zdata.row(j))% (1-Zdata.row(Zneighbor(j)-1)));
       
                   
        }
    }
    
    // storage
    return Rcpp::List::create(Rcpp::Named("tau_t") = output_t,Rcpp::Named("tau_c") = output_c );



}


