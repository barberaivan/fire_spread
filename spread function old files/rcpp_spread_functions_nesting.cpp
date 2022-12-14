#include <Rcpp.h>
using namespace Rcpp;


// parece que para que la función esté disponible para Cpp, no hay que ponerle
// Rcpp:export. Eso hace que no sea usable en R, pero quizás con Rcpp::interfaces
// sí se pueda.


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
    double linpred = coef(0); // linear predictor, initialize as intercept. (the design matrix lacks the intercept column.)
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


// Adjacent function
// [[Rcpp::export]]
IntegerMatrix adjacent_cpp(IntegerVector cells, IntegerVector n_rowcol,
                           IntegerMatrix moves) {
  
  // get row and col from cell id
  IntegerMatrix row_col = cell_to_rowcol_cpp(as<NumericVector>(cells), 
                                             as<NumericVector>(n_rowcol));
  // cast as numeric because that function needs numeric inputs.
  
  // CAREFUL: we still use R indexing (starts at 1); that will have to change.
  
  // Neighbours cells
  IntegerMatrix neigh_cells(cells.length(), 8); // [burning_cells, neighbours]
  
  for(int c = 0; c < cells.length(); c++) {
    // get neighbours row_col for cell c
    IntegerMatrix neigh_rc(2, 8); // [c(row, col), neighbours]
    
    // (loop over 8 neighbours)
    for(int n = 0; n < 8; n++) {
      neigh_rc(_, n) = row_col(_, c) + moves(_, n); // neighbours row-column pairs
    }      
    
    //Rcout << "neigh_rc for cell " << c << ": " << neigh_rc << "\n";  
    
    // Define valid neighbours values (not out of range)
    IntegerVector valid_id(8);

    int lower = 1; // change by 0 later, it's the lowest allowed row or column index.
    for(int n = 0; n < 8; n++) {
      if(neigh_rc(0, n) >= lower & neigh_rc(0, n) <= n_rowcol(0) &  // check rows
         neigh_rc(1, n) >= lower & neigh_rc(1, n) <= n_rowcol(1)) { // check columns
        valid_id(n) = 1;
      }
    }
    
    //Rcout << "valid_id filled for cell " << c << ": \n" << valid_id << "\n";  
    
    
    // Get raw neighbours row-cols (includes invalid positions)
    //IntegerVector neigh_cells_raw = rowcol_to_cell_cpp(neigh_rc, 
    //                                                   as<NumericVector>(n_rowcol));
    
    neigh_cells(c, _) = rowcol_to_cell_cpp(neigh_rc,
                                           as<NumericVector>(n_rowcol));
    
    // Make NA the cell ID for invalid neighbours
    for(int n = 0; n < 8; n++) {
      if(valid_id(n) == 0) 
        neigh_cells(c, n) = NA_INTEGER;
    }
  }
   
  return neigh_cells;
}


/*
# How to use: 
neighbours_matrix <- adjacent_vector(cells = burning_cells, n_rowcol = n_rowcol) 
  
# Function
  adjacent_vector <- function(cells, n_rowcol) {
    
# (data to test)
# cells = rep(burning_cells, 3)
# cells = burning_cells
# n_rowcol = rc
    
# get row and col from cell id
    row_col <- cell_to_rowcol(cells, n_rowcol)
      
# neighbours row_col 
      neigh_rc <- array(NA, dim = c(2, 8, length(cells)))
        for(i in 1:length(cells)) {
          neigh_rc[, , i] <- row_col[, i] + moves # neighbours row-column pairs
        }                    
        
# Get position of values out of range
        valid_id <- matrix(0, length(cells), 8)
          for(c in 1:length(cells)) {
            for(i in 1:8) {
              if(neigh_rc[1, i, c] > 0 & neigh_rc[1, i, c] <= n_rowcol[1] &
                 neigh_rc[2, i, c] > 0 & neigh_rc[2, i, c] <= n_rowcol[2]) {
                valid_id[c, i] <- 1
              }
            }
          }
          
# get cell id
          neigh_cell <- matrix(NA, length(cells), 8)
            for(c in 1:length(cells)) {
              for(i in 1:8) {
                if(valid_id[c, i] == 1) {
                  tmp <- matrix(neigh_rc[, i, c], 2, 1)
                  neigh_cell[c, i] <- rowcol_to_cell(tmp, n_rowcol)
                } else {
                  neigh_cell[c, i] <- NA
                }
              }
            }
            
            return(neigh_cell)
  }
*/