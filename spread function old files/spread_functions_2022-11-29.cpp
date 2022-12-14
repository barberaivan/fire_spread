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
 *  adjacent_cpp
 *  
 * Functions to actually spread fire:
 *  spread_around_cpp
 *  simulate_fire_cpp  
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
    // row
    rowcols(0, i) = ceil(cells[i] / n_rowcol[1]);
        // (This division needs numeric inputs to return a real value)
    // column
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
    rowcols(1, i) = n_rowcol[1] - 
                    ( (rowcols(0, i) + 1) * n_rowcol[1] - (cells[i] + 1) ) - 
                    1;
    // cells also start at 0!
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
//'   focal_row -1,
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
    neigh_cells(c, _) = rowcol_to_cell_cpp0(neigh_rc,
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
//'   Wind and elevation columns must be the last 2.
//' @param NumericVector distances: distances (m) between burning and target cells,
//'   in the same order as positions. Used to compute the elevation effect.
//'   This vector depends on the neighbourhood design and on the pixel scale.
//'   If unchanged, it's always the same.
//' @param NumericVector angles: angles (radians) between burning and target cells,
//'   in the same order as positions. Used to compute the wind effect. This
//'   vector is always the same, unless the neighborhood is changed.
//' @param double upper_limit: upper limit for spread probability (setting to 
//'   1 makes absurdly large fires). 

/*
 * Originally this functions modified the neighbours' data matrix, in the wind
 * and elevation columns, to compute angles or whatever was needed. That was
 * a bit problematic because then, in R the matrix remained modified, and 
 * consecutive runs of the function were not the same.
 * This was resolved by creating new objects slope_term and wind_term, which
 * are added to the linear predictor after the loop adds the other terms.
 * This is why the wind and elevation columns in the data must be the last ones.
 */

// [[Rcpp::export]]
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
  
  // Loop over neighbours: compute p and sample fire spread.
  for(int i = 0; i < data_neighbours.nrow(); i++) {
    
    // recompute wind and elevation columns (will be wind and slope effects)
    
    // slope (from elevation)
    double slope_term = sin(atan(
      (data_neighbours(i, elev_column) - data_burning(elev_column)) / 
        distances(positions(i))
    ));
    
    //Rcout << data_neighbours(i, elev_column) << "\n";
    
    // wind term
    double wind_term = cos(angles(positions(i)) - data_burning(wind_column));
    
    // compute probability
    double linpred = coef(0); // linear predictor, initialize as intercept. 
    // (the design matrix lacks the intercept column.)
    // All the effects besides wind and slope: 
    for(int k = 0; k < (data_neighbours.ncol() - 2); k++) {
      linpred += coef(k+1) * data_neighbours(i, k);
    }
    // (coef is lagged ahead wrt the data_neighbours because it has the intercept,
    // which is not present in the data.)
    // Add wind and slope terms
    linpred += slope_term * coef[elev_column + 1] + 
               wind_term * coef[wind_column + 1];
    
    double prob = upper_limit / (1 + exp(-linpred));
    
    // ---------------------------------------------------------------
    
    // Simulate the spread
    burn(i) = R::rbinom(1.0, prob);
  }
  
  return burn; 
  // In this way it returns a vector containing the cell id of burning neighbours.
  // return cell_neighbours[burn == 1];
}




// The same function but returning prob ----

// [[Rcpp::export]]
NumericVector spread_around_prob_cpp(NumericVector data_burning, 
                                     NumericMatrix data_neighbours,  
                                     NumericVector coef,
                                     IntegerVector positions,
                                     int wind_column,
                                     int elev_column,
                                     NumericVector distances,
                                     NumericVector angles,
                                     double upper_limit = 1.0) {

  // burned or not vector to fill
  NumericVector prob(data_neighbours.nrow());    
  
  // Loop over neighbours: compute p and sample fire spread.
  for(int i = 0; i < data_neighbours.nrow(); i++) {
    
    // recompute wind and elevation columns (will be wind and slope effects)
    
    // slope (from elevation)
    double slope_term = sin(atan(
      (data_neighbours(i, elev_column) - data_burning(elev_column)) / 
        distances(positions(i))
    ));
    
    // wind term
    double wind_term = cos(angles(positions(i)) - data_burning(wind_column));
    
    // compute probability
    double linpred = coef(0); // linear predictor, initialize as intercept. 
    // (the design matrix lacks the intercept column.)
    // All the effects besides wind and slope: 
    for(int k = 0; k < (data_neighbours.ncol() - 2); k++) {
      linpred += coef(k+1) * data_neighbours(i, k);
    }
    // (coef is lagged ahead wrt the data_neighbours because it has the intercept,
    // which is not present in the data.)
    // Add wind and slope terms
    linpred += slope_term * coef[elev_column + 1] + 
      wind_term * coef[wind_column + 1];
    
    prob(i) = upper_limit / (1 + exp(-linpred));
  }
  
  return prob; // Used to check whether it matches the R function
}






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
//' @param IntegerVector ignition_cell: id for the cell(s) where the fire begun.
//' @param IntegerVector burnable: vector indicating if each pixel is burnable (1)
//'   or not (0).
//' @param NumericVector n_rowcol: number or row and columns of the landscape.
//' @param NumericVector coef: parameters in logistic regression to compute the
//'   spread probability as a function of covariates.
//' @param IntegerVector positions: position of each neighbour in 
//'   relation to the burning cell. The eight neighbours are labelled from 0 to
//'   7 beggining from the upper-left corner (by row):
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
//'   Wind and elevation columns must be the last 2.
//' @param NumericVector distances: distances (m) between burning and target cells,
//'   in the same order as positions. Used to compute the elevation effect.
//'   This vector depends on the neighbourhood design and on the pixel scale.
//'   If unchanged, it's always the same.
//' @param NumericVector angles: angles (radians) between burning and target cells,
//'   in the same order as positions. Used to compute the wind effect. This
//'   vector is always the same, unless the neighborhood is changed.
//' @param double upper_limit: upper limit for spread probability (setting to 
//'   1 makes absurdly large fires). 

// [[Rcpp::export]]
IntegerVector simulate_fire_cpp(NumericMatrix landscape, 
                                IntegerVector ignition_cells,//    
                                IntegerVector burnable,      // burn layer in {0, ..., 3}, see code below
                                IntegerVector n_rowcol,
                                NumericVector coef,
                                IntegerMatrix moves,
                                int wind_column,
                                int elev_column,
                                NumericVector distances,
                                NumericVector angles,
                                double upper_limit) {
  
  Rcout << "ignition cells " << ignition_cells << "\n";
  
  int n_row = n_rowcol[0];
  int n_col = n_rowcol[1];
  int n_cell = n_row * n_col;
  IntegerVector cell_ids = seq(0, n_cell - 1);
  
  Rcout << "cell ids " << cell_ids << "\n";
  
  
  // Create burn layer, which will be exported.
  IntegerVector burn = rep(0, n_cell);
  
  // Fill the burn vector with 3 in the non-burnable pixels
  burn[burnable == 0] = 3;
  // Set to 1 the ignition point
  burn[ignition_cells] = 2;
  
  // Get burning cells (initialized at the ignition point).
  // This object would need variable length, but we fix it at 0.6 * n_cell
  // and record which values we have to fill. In this way, there is always enough
  // space and memory assignment is not repeated.
  int n_cell_23 = ceil(n_cell / 3);
  IntegerVector burning_cells(n_cell_23);
  
  // IntegerVector burning_cells = ignition_cells;
  Rcout << "burning_cells start " << burning_cells << "\n";
  
  // Position AFTER The last burning cell in the burning_cells vector. If 0, 
  // it means there are no burning cells anymore.
  int burning_last = ignition_cells.length();
  burning_cells[Range(0, burning_last - 1)] = ignition_cells;
  
  /*
   * Burn code: 
   * 0 = not burned but burnable,
   * 1 = burning,
   * 2 = burned, 
   * 3 = not burnable.
   */
  
  int burning_size = ignition_cells.length();
  int j = 1;
  
  // while(burning_cells.length() > 0) {
  while(burning_last > 0) {
  // while(j < 11) {  
    
    Rcout << "burn iteration " << j << "\n";  
    IntegerVector burning_start = burning_cells[Range(0, burning_last - 1)];
    Rcout << "burning cells start " << burning_start << "\n";  
    
    // IntegerVector burn_2 = burn[burn == 2];
    // Rcout << "ya quemadas " << burn_2 << "\n";
    
    // Get burning_cells' neighbours
    IntegerMatrix neighbours = adjacent_cpp0(burning_start, //burning_cells[Range(0, burning_last - 1)], 
                                             n_rowcol, moves);
    Rcout << "neighbours matrix " << neighbours << "\n";  
    
    // Loop to spread fire from every burning_cell
    // for(int b = 0; b < burning_cells.length(); b++) {
    for(int b = 0; b < burning_last; b++) {
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
      
      // If there are burnable neighbours, spread around
      if(neigh_b.length() > 0) {
        
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
        IntegerVector temp_burn = burn[neigh_b];
        Rcout << "temp_burn " << temp_burn << "\n";  
        
      }
    } // end loop across burning cells
    

    // Update burning cells ids.
    
    // double burning_last_real;
    // burning_last_real = sum(as<NumericVector>(burn[burn == 1]));
    // burning_last = as_integer(burning_last_real);
    burning_last = sum(as<NumericVector>(burn[burn == 1]));
    
    // IntegerVector temp = cell_ids[burn == 1]; // to avoid double subsetting
    IntegerVector temp = cell_ids[burn == 1];
    burning_cells[Range(0, burning_last - 1)] = temp;
    // burning_size = sum(as<NumericVector>(burn[burn == 1]));
    // burning_size = burning_cells.length();
    
    Rcout << "burn iteration " << j << "\n";  
    Rcout << "burning cells end " << temp << "\n";  
    // Rcout << "las quemando " << las_quemando << "\n";   same as burning-cells-end
    Rcout << "burning_size " << burning_last << "\n";  
    
    burn[burning_cells[Range(0, burning_last - 1)]] = 2;
    
    j += 1;
    
    // burning_cells[Rcpp::Range(0, burning_last - 1)] = temp;
    // parece que no le gusta que la asignaci??n tenga doble subsetting
    // burning_cells = cell_ids[burn == 1];
  } // end while 
  
  return burn;
}


