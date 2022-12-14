#include <Rcpp.h>
using namespace Rcpp;

/*
 * Functions to spread fire. 
 * 
 * The fire spread model is a cellular automata that spreads towards the 8-
 * neighbouring pixels. Rasters are treated as vectors, with a cell id by pixel. 
 * We need functions to translate between the vectorized cell id to row-column 
 * and the other way around.
 * 
 * In this script all functions are created for 0-starting-indexes.
 * I still keep the 1-starting indexes functions because I think that makes 
 * the code comparable with R. The 0 ending in the function names indicates
 * that it's written in 0-starting indexes.
 * 
 * Functions to find neighbours of each target cell: 
 *  cell_to_rowcol_cpp
 *  rowcol_to_cell_cpp
 *  
 *  
 */

// -----------------------------------------------------------------------

//' @title cell_to_rowcol_cpp
//' @description Translates cell id to row and column.
//' @return IntegerMatrix(2, n_cells): matrix with row and column ids, with each
//'   column corresponding to a cell.
//' @param NumericVector cells: cell ids.
//' @param NumericVector n_rowcol: number or row and columns of the landscape.

// 1 start
// [[Rcpp::export]]
IntegerMatrix cell_to_rowcol_cpp(NumericVector cells, NumericVector n_rowcol) {
  int n_cells = cells.size();
  
  IntegerMatrix rowcols(2, n_cells); // rows positions in row 0, cols positions in row 1.
  
  for(int i = 0; i < n_cells; i++) {
    // This division needs numeric inputs to return a real value
    rowcols(0, i) = ceil(cells[i] / n_rowcol[1]);
    rowcols(1, i) = n_rowcol[1] - (rowcols(0, i) * n_rowcol[1] - cells[i]);
  }
  
  return rowcols;
}

// 0 start
// [[Rcpp::export]]
IntegerMatrix cell_to_rowcol_cpp0(NumericVector cells, NumericVector n_rowcol) {
  int n_cells = cells.size();
  
  IntegerMatrix rowcols(2, n_cells); // rows positions in row 0, cols positions in row 1.
  
  for(int i = 0; i < n_cells; i++) {
    // This division needs numeric inputs to return a real value
    rowcols(0, i) = ceil((cells[i] + 1) / n_rowcol[1]) - 1;
    rowcols(1, i) =  cells[i] - rowcols(0, i) * n_rowcol[1];
  }
  
  return rowcols;
}

// -----------------------------------------------------------------------

//' @title rowcol_to_cell_cpp
//' @description Translates row and column id to cell id.
//' @return IntegerVector(n_cells): matrix with row and column ids, with each
//'   column corresponding to a cell.
//' @param IntegerMatrix rowcol: matrix with row and column ids, with each
//'   column corresponding to a cell.
//' @param NumericVector n_rowcol: number or row and columns of the landscape.

// 1 start
// [[Rcpp::export]]
IntegerVector rowcol_to_cell_cpp(IntegerMatrix rowcol, NumericVector n_rowcol) {
  
  IntegerVector cells(rowcol.ncol());
  
  for(int i = 0; i < cells.size(); i++) {
    cells(i) = (rowcol(0, i) - 1) * n_rowcol(1) + rowcol(1, i);
  }
  return cells;
}

// 0 start
// [[Rcpp::export]]
IntegerVector rowcol_to_cell_cpp0(IntegerMatrix rowcol, NumericVector n_rowcol) {
  
  IntegerVector cells(rowcol.ncol());
  
  for(int i = 0; i < cells.size(); i++) {
    cells(i) = rowcol(0, i) * n_rowcol(1) + rowcol(1, i);
  }
  return cells;
}

// -----------------------------------------------------------------------

//' @title adjacent_cpp
//' @description Gets the cell id of the neighbours of focal (burning) cells. It
//'   takes into account that neighbours can't fall outside the landscape.
//' @return IntegerMatrix(number of focal cells, 8): focal cells in rows, 
//'   neighbours in columns. (8-pixels neighbourhood).
//' @param IntegerVector cells: focal (burning) cell ids.
//' @param NumericVector n_rowcol: number or row and columns of the landscape.
//' @param IntegerMatrix moves: matrix defining the neighbourhood. In the case
//'   of a 8-pixels neighbourhood, it's a matrix with 2 rows (row and column 
//'   values) and 8 columns. Its values are {-1, 0, 1}, so that when adding up
//'   the row-col ids of a cell and a column of moves, we get the row-col of a
//'   neighbour. They are oredered like this:
//'   1 2 3
//'   4   5
//'   6 7 8
//'   For example, to get the neighbour 1, we compute 
//'   focal_row - 1,
//'   focal_column - 1.
//'   The first column in moves is then 
//'   -1
//'   -1
//'   (rows and columns are numbered from the upper-left corner.)

// Adjacent function start 1
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
    
    // Write raw neighbours to result matrix (includes invalid ones)
    neigh_cells(c, _) = rowcol_to_cell_cpp(neigh_rc,
                                           as<NumericVector>(n_rowcol));
    
    // Define valid neighbours values (not out of range)
    IntegerVector valid_id(8);
    
    int lower = 1; // change by 0 later, it's the lowest allowed row or column index.
    for(int n = 0; n < 8; n++) {
      if(neigh_rc(0, n) >= lower & neigh_rc(0, n) <= n_rowcol(0) &  // check rows
         neigh_rc(1, n) >= lower & neigh_rc(1, n) <= n_rowcol(1)) { // check columns
        valid_id(n) = 1;
      }
    }

    // Make NA the cell ID for invalid neighbours
    for(int n = 0; n < 8; n++) {
      if(valid_id(n) == 0) 
        neigh_cells(c, n) = NA_INTEGER;
    }
  }
  
  return neigh_cells;
}

// Adjacent function start 0
// [[Rcpp::export]]
IntegerMatrix adjacent_cpp0(IntegerVector cells, IntegerVector n_rowcol,
                            IntegerMatrix moves) {
  
  // get row and col from cell id
  IntegerMatrix row_col = cell_to_rowcol_cpp0(as<NumericVector>(cells), 
                                              as<NumericVector>(n_rowcol));
  // cast as numeric because that function needs numeric inputs.

  // Neighbours cells
  IntegerMatrix neigh_cells(cells.length(), 8); // [burning_cells, neighbours]
  
  for(int c = 0; c < cells.length(); c++) {
    // get neighbours row_col for cell c
    IntegerMatrix neigh_rc(2, 8); // [c(row, col), neighbours]
    
    // (loop over 8 neighbours)
    for(int n = 0; n < 8; n++) {
      neigh_rc(_, n) = row_col(_, c) + moves(_, n); // neighbours row-column pairs
    }      
    
    // Write raw neighbours to result matrix (includes invalid ones)
    neigh_cells(c, _) = rowcol_to_cell_cpp(neigh_rc,
                                           as<NumericVector>(n_rowcol));
    
    // Define valid neighbours values (not out of range)
    IntegerVector valid_id(8);
    
    int lower = 0; // it's the lowest allowed row or column index (c++ indexing)
    for(int n = 0; n < 8; n++) {
      if(neigh_rc(0, n) >= lower & neigh_rc(0, n) < n_rowcol(0) &  // check rows
         neigh_rc(1, n) >= lower & neigh_rc(1, n) < n_rowcol(1)) { // check columns
         valid_id(n) = 1; // label as 1 the valid ones
      }
    }
    
    // Make NA the cell ID for invalid neighbours
    for(int n = 0; n < 8; n++) {
      if(valid_id(n) == 0) 
        neigh_cells(c, n) = NA_INTEGER;
    }
  }
  
  return neigh_cells;
}


// -----------------------------------------------------------------------

//' @title spread_around_cpp
//' @description Spreads fire from one burning cell to its neighbours. Most of 
//'   the work implies computing the burn probability given the covariates. This
//'   is the function that should be edited when alternative spread functions 
//'   are considered. After the spread probability is computed, the only 
//'   remaining step is to simulate a Bernoulli process in each neighbour.
//' @return IntegerVector(number of target neighbours): binary vector indicating 
//'   burned (1) or not (0) for every analyzed neighbour. Note that invalid
//'   or already burned neighbours should not be included, but that subsetting
//'   is achieved by another function.

//' @param NumericVector data_burning: environmental data from burning cell
//'   (the most important here are wind direction and elevation).
//' @param NumericMatrix data_neighbours: environmental data from target 
//'   neighbours with a column by landscape layer.
//' @param NumericVector coef: parameters in logistic regression to compute the
//'   spread probability as a function of covariates.
//' @param IntegerVector positions: relative position of each neighbour in 
//'   relation to the burning cell. The eight neighbours are labelled from 0 to
//'   7 beggining from the upper-left one (by row):
//'   0 1 2
//'   3   4
//'   5 6 7.
//'   This is necessary to compute the elevation and wind effects, as they 
//'   depend on the angle and distance between burning and target pixels. 
//'   (Sometimes not all neighbours are burnable, so this vector will not
//'   always be complete.)
//' @param int wind_column: column in the data (landscape) with wind value.
//' @param int elev_column: column in the data (landscape) with elevation value.
//' @param NumericVector distances: distances (m) between burning and target cells,
//'   in the same order as positions. Used to compute the elevation effect.
//'   This vector depends on the neighbourhood design and on the pixel scale.
//'   If unchanged, it's always the same.
//' @param NumericVector angles: angles (°) between burning and target cells,
//'   in the same order as positions. Used to compute the wind effect. This
//'   vector is always the same, unless the neighborhood is changed.
//' @param double upper_limit: upper limit for spread probability (setting to 
//'   1 makes absurdly large fires). 
 
IntegerVector spread_around_cpp(NumericVector data_burning, 
                                NumericMatrix data_neighbours,  
                                NumericVector coef,
                                IntegerVector positions,
                                int wind_column,
                                int elev_column,
                                NumericVector distances,
                                NumericVector angles,
                                double upper_limit = 1.0) {
  
  // burned or not vector to fill
  IntegerVector burn(data_neighbours.nrow());    
  //NumericVector prob(data_neighbours.nrow());    
  
  // Loop over neighbours: compute p and sample fire spread.
  for(int i = 0; i < data_neighbours.nrow(); i++) {
    
    // recompute wind and elevation columns (will be wind and slope effects)
    
    // elevation (slope)
    data_neighbours(i, elev_column) = sin(atan(
      (data_neighbours(i, elev_column) - data_burning(elev_column)) / 
        distances(positions(i))
    ));
    
    //Rcout << data_neighbours(i, elev_column) << "\n";
    
    // wind
    data_neighbours(i, wind_column) = cos(angles(positions(i)) - data_burning(wind_column));
    
    // compute probability
    double linpred = coef(0); // linear predictor, initialize as intercept. 
                              // (the design matrix lacks the intercept column.)
    for(int k = 0; k < data_neighbours.ncol(); k++) {
      linpred += coef(k+1) * data_neighbours(i, k);
    }
    
    double prob = upper_limit / (1 + exp(-linpred));
    
    // ---------------------------------------------------------------
    
    // Simulate the spread
    burn(i) = R::rbinom(1.0, prob);
    /* TO DO: check whether a runif strategy is faster */
  }
  
  //return prob; // Used to check whether it matches the R function
  return burn;
  
  // In this way it returns a vector containing the cell id of burning neighbours.
  // return cell_neighbours[burn == 1];
}




// -----------------------------------------------------------------------

// NO USO ESTA FUNCIÓN


//' @title spread_update_cpp
//' @description Moves the spread process one step forward: it takes a landscape
//'   and the burn vector and simulates one spread stage from the burning cells.
//' @return IntegerVector(n_row * n_col): updated burn layer, coded as 
//'   0 = not burned but burnable,
//'   1 = burning,
//'   2 = burned, 
//'   3 = not burnable.

//' @param NumericMatrix landscape: environmental data from the whole landscape.
//' @param NumericMatrix burn: vector with burn-state of all cells, resulting
//'   from a previoius spread fire or from a fire beginning.
//' @param NumericVector n_rowcol: number or row and columns of the landscape.
//' @param NumericVector coef: parameters in logistic regression to compute the
//'   spread probability as a function of covariates.
//' @param IntegerVector positions: relative position of each neighbour in 
//'   relation to the burning cell. The eight neighbours are labelled from 0 to
//'   7 beggining from the upper-left one (by row):
//'   0 1 2
//'   3   4
//'   5 6 7.
//'   This is necessary to compute the elevation and wind effects, as they 
//'   depend on the angle and distance between burning and target pixels. 
//'   (Sometimes not all neighbours are burnable, so this vector will not
//'   always be complete.)
//' @param int wind_column: column in the data (landscape) with wind value.
//' @param int elev_column: column in the data (landscape) with elevation value.
//' @param NumericVector distances: distances (m) between burning and target cells,
//'   in the same order as positions. Used to compute the elevation effect.
//'   This vector depends on the neighbourhood design and on the pixel scale.
//'   If unchanged, it's always the same.
//' @param NumericVector angles: angles (°) between burning and target cells,
//'   in the same order as positions. Used to compute the wind effect. This
//'   vector is always the same, unless the neighborhood is changed.
//' @param double upper_limit: upper limit for spread probability (setting to 
//'   1 makes absurdly large fires). 

// [[Rcpp:export]]
// IntegerVector spread_update_cpp(NumericMatrix landscape, 
//                                 IntegerVector burned,    // burn layer in {0, ..., 3}
//                                 IntegerVector n_rowcol,
//                                 NumericVector coef,
//                                 IntegerMatrix moves,
//                                 int wind_column,
//                                 int elev_column,
//                                 NumericVector distances,
//                                 NumericVector angles,
//                                 double upper_limit) {
//   
//   int n_row = n_rowcol[0];
//   int n_col = n_rowcol[1];
//   
//   IntegerVector burning_cells = burned[burned == 1];    
//   
//   /*
//    * Burn code: 
//    * 0 = not burned but burnable,
//    * 1 = burning,
//    * 2 = burned, 
//    * 3 = not burnable.
//    */
//   
//   // Get burning_cells' neighbours
//   IntegerMatrix neighbours = adjacent_cpp0(burning_cells, n_rowcol, moves);
//   
//   // Loop to spread fire from every burning_cell
//   for(int b = 0; b < burning_cells.length(); b++) {
//     // Get burning_cells' data
//     NumericVector data_burning = landscape(burning_cells[b], _);
//     
//     // Get its valid neighbours (not NA) and record their relative positions
//     // (labelled as follows: 
//     // 0 1 2 
//     // 3   4
//     // 5 6 7)
//     // This is necessary to compute the elevation and wind effects, as they 
//     // depend on the angle and distance between burning and target pixels.
//     IntegerVector neigh_b = neighbours(b, _);
//     IntegerVector positions = seq(0, 7);
//     // Make NA if the pixel is not burnable (not 0)
//     for(int n = 0; n < 8; n++) {
//       if((neigh_b[n] == NA_INTEGER) | (burned[neigh_b[n]] != 0)) {
//         neigh_b[n] = NA_INTEGER;
//         positions[n] = NA_INTEGER;
//       }
//     }
//     // remove NAs
//     neigh_b = na_omit(neigh_b);
//     positions = na_omit(positions);
//     
//     // Get its valid neighbours' data 
//     NumericMatrix data_neighbours(neigh_b.length(), landscape.ncol());
//     for(int n = 0; n < neigh_b.length(); n++) {
//       data_neighbours(n, _) = landscape(neigh_b[n], _);
//     }
//     
//     // Spread fire towards neighbours and write to burned vector
//     burned[neigh_b] = spread_around_cpp(data_burning, 
//                                         data_neighbours,  
//                                         cell_neighbours = neigh_b,
//                                         coef,
//                                         positions, // columns of valid neighbours
//                                         wind_column,
//                                         elev_column,
//                                         distances,
//                                         angles,
//                                         upper_limit);
//     
//   }
//   
//   // Update step: which was previously burning (1) is now burned (2).
//   burned[burning_cells] = 2;
//   return burned; // replaces the previous burned vector, which was provided
//                  // as data.
// } 
// 
// // lo que sigue es testear esto y luego agregar el while.



// -----------------------------------------------------------------------

// function to simulate a fire spread given the landscape, 
// model coefficients and ignition points.

//' @title simulate_fire_cpp
//' @description function to simulate a fire spread given the landscape, 
//'   model coefficients and ignition points.
//' @return IntegerVector(n_row * n_col): updated burn layer, coded as 
//'   0 = not burned but burnable,
//'   1 = burning (only occurs before the function runs),
//'   2 = burned, 
//'   3 = not burnable.

//' @param NumericMatrix landscape: environmental data from the whole landscape.
//' @param NumericMatrix burn: vector with burn-state of all cells. The ignition
//'   points are coded as ones. This function will update this layer and return
//'   a burned one, where there should be no ones.
//' @param NumericVector n_rowcol: number or row and columns of the landscape.
//' @param NumericVector coef: parameters in logistic regression to compute the
//'   spread probability as a function of covariates.
//' @param IntegerVector positions: relative position of each neighbour in 
//'   relation to the burning cell. The eight neighbours are labelled from 0 to
//'   7 beggining from the upper-left one (by row):
//'   0 1 2
//'   3   4
//'   5 6 7.
//'   This is necessary to compute the elevation and wind effects, as they 
//'   depend on the angle and distance between burning and target pixels. 
//'   (Sometimes not all neighbours are burnable, so this vector will not
//'   always be complete.)
//' @param IntegerMatrix moves: see adjacent_cpp0 documentation.
//' @param int wind_column: column in the data (landscape) with wind value.
//' @param int elev_column: column in the data (landscape) with elevation value.
//' @param NumericVector distances: distances (m) between burning and target cells,
//'   in the same order as positions. Used to compute the elevation effect.
//'   This vector depends on the neighbourhood design and on the pixel scale.
//'   If unchanged, it's always the same.
//' @param NumericVector angles: angles (°) between burning and target cells,
//'   in the same order as positions. Used to compute the wind effect. This
//'   vector is always the same, unless the neighborhood is changed.
//' @param double upper_limit: upper limit for spread probability (setting to 
//'   1 makes absurdly large fires). 

IntegerVector simulate_fire_cpp(NumericMatrix landscape, 
                                IntegerVector burn,    // burn layer in {0, ..., 3}, with 1s in the ignition points
                                IntegerVector n_rowcol,
                                NumericVector coef,
                                IntegerMatrix moves,
                                int wind_column,
                                int elev_column,
                                NumericVector distances,
                                NumericVector angles,
                                double upper_limit) {
  
  int n_row = n_rowcol[0];
  int n_col = n_rowcol[1];
  int n_cell = n_row * n_cell;
  IntegerVector cell_ids = seq(0, n_cell - 1);
  
  // Get burning cells (initialized at the ignition point)
  IntegerVector burning_cells = cell_ids[burn == 1];

  
  /*
   * Burn code: 
   * 0 = not burned but burnable,
   * 1 = burning,
   * 2 = burned, 
   * 3 = not burnable.
   */
  
  while(burning_cells.length() > 0) {
    
    // Get burning_cells' neighbours
    IntegerMatrix neighbours = adjacent_cpp0(burning_cells, n_rowcol, moves);
    
    // Loop to spread fire from every burning_cell
    for(int b = 0; b < burning_cells.length(); b++) {
      // Get burning_cells' data
      NumericVector data_burning = landscape(burning_cells[b], _);
      
      // Get its valid neighbours (not NA) and record their relative positions
      // (labelled as follows: 
      // 0 1 2 
      // 3   4
      // 5 6 7)
      // This is necessary to compute the elevation and wind effects, as they 
      // depend on the angle and distance between burning and target pixels.
      IntegerVector neigh_b = neighbours(b, _);
      IntegerVector positions = seq(0, 7);
      // Make NA if the pixel is not burnable (not 0) or is not a valid neighbour
      for(int n = 0; n < 8; n++) {
        if(R_IsNA(neigh_b[n]) | (burn[neigh_b[n]] != 0)) {
          neigh_b[n] = NA_INTEGER;
          positions[n] = NA_INTEGER;
        }
      }
      // remove NAs
      neigh_b = na_omit(neigh_b);
      positions = na_omit(positions);
      
      // Get its valid neighbours' data 
      NumericMatrix data_neighbours(neigh_b.length(), landscape.ncol());
      for(int n = 0; n < neigh_b.length(); n++) {
        data_neighbours(n, _) = landscape(neigh_b[n], _);
      }
      
      // Spread fire towards neighbours and write to burned vector
      burn[neigh_b] = spread_around_cpp(data_burning, 
                                        data_neighbours,  
                                        coef,
                                        positions, // columns of valid neighbours
                                        wind_column,
                                        elev_column,
                                        distances,
                                        angles,
                                        upper_limit);
    
    } // end loop across burning cells
  
    // Turn burning into burned (2)
    burn[burning_cells] = 2;
    // Update burning cells ids.
    burning_cells = cell_ids[burn == 1];
  } // end while 
  
  return burn;
}


// ----------------------------------------------------------------------


// sigue testear en R spread_around y simulate_fire. 
// Me suena que ya había probado la de spread around...
