#include <RcppArmadillo.h>
#include <random>


// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// [[Rcpp::export()]]

Rcpp::List calculate_tau(arma::mat Y, arma::mat Zdata){

    // Define variables
    int Ny = Y.n_rows;
    int Nz = Zdata.n_rows;

    // Final outoput matrix
    arma::mat output(Ny, Nz);
    
    for ( int i=0; i<Ny ; i++){

        for (int j=0; j<Nz; j++){

          
           output(i,j) = arma::sum(Y.row(i) % Zdata.row(j)) / arma::sum (Zdata.row(j)) - arma::sum(Y.row(i) % (1-Zdata.row(j)))/ arma::sum (1-Zdata.row(j));
        
                   
        }
    }
    
    // storage
    return Rcpp::List::create(Rcpp::Named("tau") = output);



}


