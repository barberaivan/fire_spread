#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix cell_to_rowcol_cpp(NumericVector cell, NumericVector n_rowcol) {
  int n_cells = cell.size();

    IntegerMatrix rowcols(2, n_cells); // rows positions in row 0, cols positions in row 1.
  
  for(int i = 0; i < n_cells; i++) {
    // This division needs numeric inputs to return a real value
    rowcols(0, i) = ceil(cell[i] / n_rowcol[1]);
    rowcols(1, i) = n_rowcol[1] - (rowcols(0, i) * n_rowcol[1] - cell[i]);
  }
  
  return rowcols;
}

// [[Rcpp::export]]
IntegerVector rowcol_to_cell_cpp(IntegerMatrix rowcol, NumericVector n_rowcol) {
  
  IntegerVector cells(rowcol.ncol());
  
  for(int i = 0; i < cells.size(); i++) {
    cells(i) = (rowcol(0, i) - 1) * n_rowcol(1) + rowcol(1, i);
  }
  return cells;
}


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
  
  // burned or not vector to fill
  NumericVector burn(data_neighbours.nrow());    
  NumericVector prob(data_neighbours.nrow());    
  
  // Loop over neighbours: compute p and sample fire spread.
  for(int i = 0; i < data_neighbours.nrow(); i++) {
    
    // recompute wind and elevation columns (will be wind and slope effects)
    
    // elevation (slope)
    data_neighbours(i, elev_column) = sin(atan(
      (data_neighbours(i, elev_column) - data_burning(elev_column)) / 
        distance(columns(i))
    ));
    
    //Rcout << data_neighbours(i, elev_column) << "\n";
    
    // wind
    data_neighbours(i, wind_column) = cos(angles(columns(i)) - data_burning(wind_column));
    
    // compute probability
    double linpred = coef(0); // linear predictor, initialize as intercept. 
                              // (the design matrix lacks the intercept column.)
    for(int k = 0; k < data_neighbours.ncol(); k++) {
      linpred += coef(k+1) * data_neighbours(i, k);
    }
    
    //double prob = 1 / (1 + exp(-linpred)) * upper_limit;
    prob(i) = 1 / (1 + exp(-linpred)) * upper_limit;
    // sample the burn
    burn(i) = R::rbinom(1.0, prob(i));
    /*
     * TO DO: check whether a runif strategy works better.
     */
  }
  
  return burn;
  //return prob; // Used to check whether it matches the R function
}

/*
 * Adjacent requires the use of internal cpp functions, so I continue in another 
 * script.
 */
