#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector spread_prob_cpp(NumericVector data_burning, 
                                  NumericMatrix data_neighbours,  
                                  NumericVector coef,
                                  IntegerVector columns,
                                  int wind_column,
                                  int elev_column,
                                  NumericVector distance,
                                  NumericVector angles,
                                  double upper_limit = 1.0) {
      
      // probability vector to fill
      NumericVector linpred(data_neighbours.nrow()); // linear predictor
      NumericVector probs(data_neighbours.nrow());   // spread probability
      NumericVector burn(data_neighbours.nrow());    // burned or not
      
      // recompute wind and elevation columns (will be wind and slope effects)
      for(int i = 0; i < data_neighbours.nrow(); i++) {
        
        // elevation (slope)
        data_neighbours(i, elev_column) = sin(atan(
          (data_neighbours(i, elev_column) - data_burning(elev_column)) / 
            distance(columns(i))
        ));
        
        // wind
        data_neighbours(i, wind_column) = cos(angles(columns(i)) - data_burning(wind_column));
        
        // compute probability
        linpred(i) = coef(0); // linear predictor, initialize as intercept. (the design matrix lacks the intercept column.)
        for(int k = 0; k < data_neighbours.ncol(); k++) {
          linpred(i) += coef(k+1) * data_neighbours(i, k);
        }
        
        probs(i) = 1 / (1 + exp(-linpred(i))) * upper_limit;
        // sample the burn
        burn(i) = R::rbinom(1.0, probs(i));
        
      }
      
      //return (probs);
      return burn;
    }
  
  