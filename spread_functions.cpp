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
 * I still keep the 1-starting indexes functions because that makes
 * the code comparable with R and prevents errors at the first try.
 * The 0 ending in the function names indicates that it's written with
 * 0-starting indexes.
 *
 * Functions to find neighbours of each target cell:
 *  cell_to_rowcol_cpp
 *  rowcol_to_cell_cpp
 *  adjacent_cpp
 *  adjacent_cpp0_2 (used this one because of a bug)
 *
 * Functions to actually spread fire:
 *  spread_around_cpp
 *    (spreads fire towards the neighbours of one burning cell. 
 *    IMPORTANT: this function holds the spread model itself, so it has to be 
 *    edited if the model is changed.)
 *  simulate_fire_cpp
 *    (basically, a while loop running spread_around_cpp as long as there
 *    are burning cells)
 *
 *
 *  NOTES: Rcpp is a package which links C++ with R. It has a few classes and
 *  functions useful to write more compact code, for example, classes as
 *  NumericMatric, IntegerVector, etc.
 *
 */

// Constants ---------------------------------------------------------------

// Elevation data to standardize predictor
double elevation_mean = 1163.3;
double elevation_sd = 399.5;

// Moves matrix to find neighbours (queen directions) of a pixel, used in
// adjacent_ functions.
IntegerVector moves_long = {-1,-1,-1,  0,0,  1,1,1,  // rows
                            -1, 0, 1, -1,1, -1,0,1}; // cols
IntegerMatrix moves_t(8, 2, moves_long.begin());
IntegerMatrix moves = transpose(moves_t);
/* In the case
 * of a 8-pixels neighbourhood, it's a matrix with 2 rows (row and column
 * values) and 8 columns. Its values are {-1, 0, 1}, so that when adding up
 * the row-col ids of a cell and a column of moves, we get the row-col of a
 * neighbour. They are oredered like this:
 * 1 2 3
 * 4   5
 * 6 7 8
 * For example, to get the neighbour 1, we compute
 * focal_row - 1,
 * focal_column - 1.
 * The first column in moves is then
 * -1
 * -1
 * (rows and columns are numbered from the upper-left corner.)
 */

// Angles between cells to compute wind effect. As the wind direction is 
// the direction from which the wind comes, these angles must represent where the
// fire would come from if from the neighbours we look at the central pixel.
NumericVector angles_raw = {
  135, 180, 225,
  90,       270,
  45,   0,  315
};
NumericVector angles = angles_raw * M_PI / 180; // in radians!


// -----------------------------------------------------------------------
  
//' @title cell_to_rowcol_cpp
//' @description Translates cell id to row and column.
//' @return IntegerMatrix(2, n_cells): matrix with row and column ids, with each
//'   column corresponding to a cell.
//' @param NumericVector cells: cell ids.
//' @param NumericVector n_rowcol: number of row and columns of the landscape.

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

// Adjacent function start 1
// [[Rcpp::export]]
IntegerMatrix adjacent_cpp(IntegerVector cells, IntegerVector n_rowcol) {

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
IntegerMatrix adjacent_cpp0(IntegerVector cells, IntegerVector n_rowcol) {

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


// Adjacent function start 0
// In this case, returns a vector, used to compute neighbours from only one 
// burning cell at the time.

// [[Rcpp::export]]
IntegerVector adjacent_vec_cpp0(IntegerVector cells, IntegerVector n_rowcol) {
  
  // get row and col from cell id
  IntegerMatrix row_col = cell_to_rowcol_cpp0(as<NumericVector>(cells),
                                              as<NumericVector>(n_rowcol));
  // cast as numeric because that function needs numeric inputs.
  
  // Neighbours cells
  IntegerMatrix neigh_cells(cells.length(), 8); // [burning_cells, neighbours]
  
  // get neighbours row_col for cell c
  IntegerMatrix neigh_rc(2, 8); // [c(row, col), neighbours]
  
  // (loop over 8 neighbours)
  for(int n = 0; n < 8; n++) {
    neigh_rc(_, n) = row_col + moves(_, n); // neighbours row-column pairs
  }
  
  // Write raw neighbours to result matrix (includes invalid ones)
  neigh_cells(0, _) = rowcol_to_cell_cpp0(neigh_rc,
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
      neigh_cells(0, n) = NA_INTEGER;
  }
  
  // Turn 1-row matrix into vector
  IntegerVector result(8);
  for(int i = 0; i < 8; i++) result[i] = neigh_cells(0, i);
  
  return result;
}



// Adjacent function start 0, using the full binary vector <burning> instead
// of burning_cells. It is used to avoid memory allocation problems.

//' @title adjacent_cpp0_2
//' @description Gets the cell id of the neighbours of focal (burning) cells. It
//'   takes into account that neighbours can't fall outside the landscape.
//' @return IntegerMatrix(number of focal cells, 9): focal cells in rows,
//'   neighbours in columns. (8-pixels neighbourhood). The first column contains
//'   the ids of the burning cells.
//' @param IntegerVector burning: binary vector as long as the landscape with
//'   1s in the cells that are currently burning.

// [[Rcpp::export]]
IntegerMatrix adjacent_cpp0_2(IntegerVector burning, IntegerVector n_rowcol) {

  int n_cells = n_rowcol[0] * n_rowcol[1];
  IntegerVector cell_ids = seq(0, n_cells - 1);
  int burning_length = sum(burning);
  IntegerVector cells(burning_length);
  cells = cell_ids[burning == 1];

  // get row and col from cell id
  IntegerMatrix row_col = cell_to_rowcol_cpp0(as<NumericVector>(cells),
                                              as<NumericVector>(n_rowcol));
  // cast as numeric because that function needs numeric inputs.

  // Neighbours cells
  IntegerMatrix neigh_cells(cells.length(), 8 + 1); // [burning_cells, c(neighbours, burning_cell id)]

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

    // Add burning cell id
    neigh_cells(c, 8) = cells[c];
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
//'   spread probability as a function of covariates. It has one more elements
//'   than the columns of data_ because it includes the intercept.
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
//' @param double upper_limit: upper limit for spread probability (setting to
//'   1 makes absurdly large fires).
   
//' The layers in data_ are:
//'   subalpine forest, {0, 1} (shrubland goes in the intercept)
//'   wet forest,       {0, 1}
//'   dry forest,       {0, 1}
//'   fwi,              (-Inf, +Inf) (standardized anomaly at pixel level)
//'   aspect            [-1, 1] # Northwestyness
//'   wind direction    [-1, 1] (windward vs leeward)
//'   elevation,        (standardized in the code)

//' The coefficients vector has a parameter by layer plus an intercept and the
//' slope effect:
//'   [Intercept] shrubland
//'   subalpine forest,
//'   wet forest,
//'   dry forest,
//'   fwi,
//'   aspect,
//'   wind,
//'   elevation,  (note slope comes after elevation)
//'   [slope],    (downhill or uphill, (-1, 1): 0 = flat, 1 = above, -1 = below)

/*
 * Originally this functions modified the neighbours' data matrix, in the wind
 * and elevation columns, to compute angles or whatever was needed. That was
 * a bit problematic because then, in R the matrix remained modified, and
 * consecutive runs of the function were not the same.
 * This was resolved by creating new objects elev_term, slope_term and wind_term, 
 * which are added to the linear predictor after the loop adds the other terms.
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
                                double upper_limit = 1.0) {

  // burned or not vector to fill
  IntegerVector burn(data_neighbours.nrow());

  // Loop over neighbours: compute p and sample fire spread.
  for(int i = 0; i < data_neighbours.nrow(); i++) {

    // compute elevation, slope and wind effects

    // slope (from elevation)
    double slope_term = sin(atan(
      (data_neighbours(i, elev_column) - data_burning(elev_column)) /
      distances(positions(i))
    ));

    // wind term
    double wind_term = cos(angles(positions(i)) - data_burning(wind_column));
    
    // elevation term (standardize predictor)
    double elev_term = (data_neighbours(i, elev_column) - elevation_mean) / elevation_sd;
    
    // compute probability
    double linpred = coef(0); // linear predictor, initialize as intercept.
    // (the design matrix lacks the intercept column.)
    // All the effects besides elevation, slope and slope:
    for(int k = 0; k < (data_neighbours.ncol() - 2); k++) {
      linpred += coef(k+1) * data_neighbours(i, k);
    }
    // (coef is lagged ahead wrt the data_neighbours because it has the intercept,
    // which is not present in the data.)
    // Add elevation, slope and wind effects
    linpred += wind_term * coef[wind_column + 1] +  
               elev_term * coef[elev_column + 1] + 
               slope_term * coef[elev_column + 2]; // slope coef is after the elevation one

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
                                     double upper_limit = 1.0) {

  // burned or not vector to fill
  NumericVector prob(data_neighbours.nrow());

  // Loop over neighbours: compute p and sample fire spread.
  for(int i = 0; i < data_neighbours.nrow(); i++) {

    // compute elevation, slope and wind effects
    
    // slope (from elevation)
    double slope_term = sin(atan(
      (data_neighbours(i, elev_column) - data_burning(elev_column)) /
        distances(positions(i))
    ));
    
    // wind term
    double wind_term = cos(angles(positions(i)) - data_burning(wind_column));
    
    // elevation term (standardize predictor)
    double elev_term = (data_neighbours(i, elev_column) - elevation_mean) / elevation_sd;
    
    // compute probability
    double linpred = coef(0); // linear predictor, initialize as intercept.
    // (the design matrix lacks the intercept column.)
    // All the effects besides elevation, slope and slope:
    for(int k = 0; k < (data_neighbours.ncol() - 2); k++) {
      linpred += coef(k+1) * data_neighbours(i, k);
    }
    // (coef is lagged ahead wrt the data_neighbours because it has the intercept,
    // which is not present in the data.)
    // Add elevation, slope and wind effects
    linpred += wind_term * coef[wind_column + 1] +   
               elev_term * coef[elev_column + 1] +
               slope_term * coef[elev_column + 2]; // slope coef is after the elevation one
               

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
//' @return IntegerVector(n_row * n_col): burned layer, coded as
//'   0 = not burned,
//'   1 = burned.

//' @param NumericMatrix landscape: environmental data from the whole landscape.
//'   See description in spread_around.
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
//' @param double upper_limit: upper limit for spread probability (setting to
//'   1 makes absurdly large fires).


// [[Rcpp::export]]
IntegerVector simulate_fire_cpp(NumericMatrix landscape,
                                IntegerVector ignition_cells,//
                                IntegerVector burnable,      // burn layer in {0, ..., 3}, see code below
                                IntegerVector n_rowcol,
                                NumericVector coef,
                                int wind_column,
                                int elev_column,
                                NumericVector distances,
                                double upper_limit) {
  
  int n_row = n_rowcol[0];
  int n_col = n_rowcol[1];
  int n_cell = n_row * n_col;
  
  // burned_ids will be filled with the cell ids of the burning pixels. start
  // and end integers will define the positions limits corresponding to the 
  // burning cells in every burn cycle.
  IntegerVector burned_ids = rep(NA_INTEGER, n_cell); // filled with integer_NA
  
  int start = 0;
  // end is the last non-empty position in the burned_ids vector.
  int end = ignition_cells.length() - 1; 
  // Example:
  // burned_ids = {231, 455, 342, 243, NA, NA, NA, NA};
  //               start          end.
  // if only one cell is burning, start = end.
  
  // initialize burned_ids and burning_size with ignition_cells
  burned_ids[seq(start, end)] = ignition_cells;
  int burning_size = ignition_cells.length(); // == end + 1 - start
  
  // The burned_bin vector will indicate whether each pixel is burned or burning
  // (1) or not (0). It's necessary to have this now because it will be used
  // to define burnable neighbours.
  IntegerVector burned_bin = rep(0, n_cell);
  // initialize with ignition_cells
  burned_bin[ignition_cells] = 1;
  
  // spread
  while(burning_size > 0) {
    // Loop over all the burning cells to burn their neighbours. Use end_forward
    // to update the last position in burned_ids within this loop, without 
    // compromising the loop's integrity.
    int end_forward = end;
    
    // b is going to keep the position in burned_ids that have to be evaluated 
    // in this burn cycle
    
    IntegerVector burned_ids_focal(burning_size);
    burned_ids_focal = burned_ids[seq(start, end)];
    
    // Rcout << "burned_ids in cycle:\n" << burned_ids_focal << "\n";
    
      
    for(int b = start; b <= end; b++) {
      
      // get id of focal burning cell
      IntegerVector focal_id(1);
      focal_id[0] = burned_ids[b];
      // Rcout << "focal_id:\n" << focal_id << "\n";
      
      // Get burning_cells' data
      NumericVector data_burning = landscape(focal_id[0], _);
      
      // get neighbours
      IntegerVector neigh_b = adjacent_vec_cpp0(focal_id, n_rowcol);
      // Rcout << "neigh_b:\n" << neigh_b << "\n";
      
      // Filter valid neighbours
      IntegerVector positions = seq(0, 7);

      // Make NA if the pixel is invalid
      for(int n = 0; n < 8; n++) {
        
        bool out_of_range = (neigh_b[n] == -2147483648); // it's NA for integers in cpp
        
        // First check if cell is out of range, because if true, we cant evaluate
        // this cell in the {burning, burnable, burned} vectors
        if(out_of_range) {
          neigh_b[n] = NA_INTEGER;
          positions[n] = NA_INTEGER;
        }
        // If it's in range, evaluate whether it's currently burnable
        else {
          bool burnable_cell = (burned_bin[neigh_b[n]] == 0) &
                               (burnable[neigh_b[n]] == 1);
          
          if(!burnable_cell) {
            neigh_b[n] = NA_INTEGER;
            positions[n] = NA_INTEGER;
          }
        }
      }
      
      // check neigh_b is ok
      // Rcout << "neigh_b:\n" << neigh_b << "\n";
      
      
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
        
        // Spread fire towards neighbours 
        IntegerVector burn_result(neigh_b.length());
        
        burn_result = spread_around_cpp(
          data_burning,
          data_neighbours,
          coef,
          positions, // columns of valid neighbours
          wind_column,
          elev_column,
          distances,
          upper_limit
        );
        
        // Rcout << "burn_result:\n" << burn_result << "\n";
        
        // store ids of recently burned cells and
        // set ones in burned_bin
        // (but advance end_forward first)
        for(int i = 0; i < neigh_b.length(); i++) {
          if(burn_result[i] == 1) {
            end_forward += 1;
            burned_ids[end_forward] = neigh_b[i];
            burned_bin[neigh_b[i]] = 1;
          }
        }
        
      } // end if for length(neighbours > 0)
    
    // Rcout << "burned_ids last:\n" << burned_ids[end_forward] << "\n";
      
      
    } // end loop over burning cells
    
    // update start and end
    start = end + 1;
    end = end_forward;
    burning_size = end - start + 1;
    
    
    // Rcout << "burning_size:\n" << burning_size << "\n";
    
    
  } // end while
  
  return burned_bin;
}


/*
 * simulate_fire_fool_cpp is probably a bad implementation of the function because
 * it evaluates the whole burn and burning vector in every burn cycle (i.e.,
 * iteration in the while loop).
 */

// [[Rcpp::export]]
IntegerVector simulate_fire_fool_cpp(NumericMatrix landscape,
                                IntegerVector ignition_cells,//
                                IntegerVector burnable,      // burn layer in {0, ..., 3}, see code below
                                IntegerVector n_rowcol,
                                NumericVector coef,
                                int wind_column,
                                int elev_column,
                                NumericVector distances,
                                double upper_limit) {

  int n_row = n_rowcol[0];
  int n_col = n_rowcol[1];
  int n_cell = n_row * n_col;
  IntegerVector cell_ids = seq(0, n_cell - 1);

  // Create burned layer, which will be exported.
  IntegerVector burned = rep(0, n_cell);
  // Set as 1 the ignition point
  burned[ignition_cells] = 1;

  // Create currently burning layer, which is now the same as burned.
  IntegerVector burning = Rcpp::clone(burned);

  // Keep size of burning pixels
  int burning_size = sum(burning);

  // while(burning_cells.length() > 0) {
  while(burning_size > 0) {
    // Get burning_cells' neighbours
    IntegerMatrix neighbours = adjacent_cpp0_2(burning,
                                               n_rowcol);

    // adjacent() needs a binary vector "burning", but once neighbours are obtained,
    // we turn the burning ones into 2, so we can differentiate new burns from
    // old ones.
    burning[burning == 1] = 2;

    // Loop to spread fire from every burning_cell
    // for(int b = 0; b < burning_cells.length(); b++) {
    for(int b = 0; b < burning_size; b++) {

      // Get burning_cells' data
      NumericVector data_burning = landscape(neighbours(b, 8), _);
      // neighbours[b, 8] holds the burning_cell id.

      // Get its valid neighbours (not NA) and record their relative positions
      // (labelled as follows:
      // 0 1 2
      // 3   4
      // 5 6 7)
      // This is necessary to compute the elevation and wind effects, as they
      // depend on the angle and distance between burning and target pixels.
      IntegerVector positions = seq(0, 7);
      IntegerVector neigh_b (8);
      for(int i = 0; i < 8; i++) neigh_b[i] = neighbours(b, i);

      // Make NA if the pixel is invalid
      for(int n = 0; n < 8; n++) {

        // bool out_of_range = R_IsNA(neigh_b[n]); // this did not work
        bool out_of_range = neigh_b[n] == -2147483648; // it's NA for integers in cpp

        // First check if cell is out of range, because if true, we cant evaluate
        // this cell in the {burning, burnable, burned} vectors
        if(out_of_range) {
          neigh_b[n] = NA_INTEGER;
          positions[n] = NA_INTEGER;
        }
        // If it's in range, evaluate whether it's currently burnable
        else {
          bool burnable_cell = (burning[neigh_b[n]] == 0) &
                               (burnable[neigh_b[n]] == 1) &
                               (burned[neigh_b[n]] == 0);

          if(!burnable_cell) {
            neigh_b[n] = NA_INTEGER;
            positions[n] = NA_INTEGER;
          }
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
        burning[neigh_b] = spread_around_cpp(
          data_burning,
          data_neighbours,
          coef,
          positions, // columns of valid neighbours
          wind_column,
          elev_column,
          distances,
          upper_limit
        );


      }
    } // end loop across burning cells

    burning[burning == 2] = 0; // turn off previously burned cells
    burned[burning == 1] = 1;  // set as burned the currently burning.

    burning_size = sum(burning);
  } // end while

  return burned;
}


// -----------------------------------------------------------------------

// function to simulate a fire spread given the landscape,
// model coefficients and ignition points.
// In this case, neighbours are precomputed, so we use more RAM but
// perhaps spend less time.
// (Made for benchmarking)

//' @title simulate_fire_cpp_notadj
//' @description function to simulate a fire spread given the landscape,
//'   model coefficients and ignition points.
//' @return IntegerVector(n_row * n_col): updated burn layer, coded as
//'   0 = not burned but burnable,
//'   1 = burning (only occurs before the function runs),
//'   2 = burned,
//'   3 = not burnable.

//' @param NumericMatrix landscape: environmental data from the whole landscape.
//' @param IntegerMatrix neighbours_matrix: matrix with the neighbour id (queen)
//'   for every cell. (n_cells rows, 8 columns)
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
//' @param double upper_limit: upper limit for spread probability (setting to
//'   1 makes absurdly large fires).

// [[Rcpp::export]]
IntegerVector simulate_fire_fool_cpp_notadj(
    NumericMatrix landscape,
    IntegerMatrix neighbours_matrix,
    IntegerVector ignition_cells,//
    IntegerVector burnable,      // burn layer in {0, ..., 3}, see code below
    IntegerVector n_rowcol,
    NumericVector coef,
    int wind_column,
    int elev_column,
    NumericVector distances,
    double upper_limit
  ) {

  int n_row = n_rowcol[0];
  int n_col = n_rowcol[1];
  int n_cell = n_row * n_col;
  int out_value = n_cell + 1;  // cell value for neighbours out of range
  IntegerVector cell_ids = seq(0, n_cell - 1);

  // Create burned layer, which will be exported.
  IntegerVector burned = rep(0, n_cell);
  // Set as 1 the ignition point
  burned[ignition_cells] = 1;

  // Create currently burning layer, which is now the same as burned.
  IntegerVector burning = Rcpp::clone(burned);

  // Keep size of burning pixels
  int burning_size = sum(burning);
  IntegerVector burning_cells;

  // while(burning_cells.length() > 0) {
  while(burning_size > 0) {

    // Rcout << "burning_size:\n" << burning_size << "\n";


    burning_cells = cell_ids[burning == 1]; // MISMO BUG QUE ANTES SI USO ESTO
    // Rcout << "burning_cells:\n" << burning_cells << "\n";

    // Get burning_cells' neighbours (NO ADJACENT EVALUATION)
    // IntegerMatrix neighbours = adjacent_cpp0_2(burning,
    //                                            n_rowcol, moves);

    IntegerMatrix neighbours(burning_size, 9); // add 9th column for burning_id pixel
    // neighbours = neighbours_matrix[burning_cells, _];

    /*
    for(int i = 0; i < burning_size; i++) {
      for(int j = 0; j < 8; j++) {
        // neighbours[i, j] = neighbours_matrix[burning_cells[i], j];
        neighbours(i, j) = neighbours_matrix[burning_cells[i], j];
      }
    }

    // Lo hago sin usar burning_cells
    */

    int counting = 0;
    for(int i = 0; i < n_cell; i++) {
      if(burning[i] == 1) {
        counting += 1;
        int row_index = counting - 1;

        // Rcout << "row_index:\n" << row_index << "\n";

        for(int j = 0; j < 8; j++) {
          // neighbours[i, j] = neighbours_matrix[burning_cells[i], j];
          // Rcout << "neighbours_matrix[i, j]:\n" << neighbours_matrix(i, j) << "\n";
          // neighbours_matrix[i, j]
          neighbours(row_index, j) = neighbours_matrix(i, j);
        }

        neighbours(row_index, 8) = i; // add burning cell id in the last column
      }
    }

    // Rcout << "neighbours:\n" << neighbours << "\n";
    // No le gusta este subseteo

    // adjacent() needs a binary vector "burning", but once neighbours are obtained,
    // we turn the burning ones into 2, so we can differentiate new burns from
    // old ones.
    burning[burning == 1] = 2;

    // Loop to spread fire from every burning_cell
    // for(int b = 0; b < burning_cells.length(); b++) {
    for(int b = 0; b < burning_size; b++) {

      // Get burning_cells' data
      NumericVector data_burning = landscape(neighbours(b, 8), _);
      // neighbours[b, 8] holds the burning_cell id.

      // Get its valid neighbours (not NA) and record their relative positions
      // (labelled as follows:
      // 0 1 2
      // 3   4
      // 5 6 7)
      // This is necessary to compute the elevation and wind effects, as they
      // depend on the angle and distance between burning and target pixels.
      IntegerVector positions = seq(0, 7);
      IntegerVector neigh_b (8);
      for(int i = 0; i < 8; i++) neigh_b[i] = neighbours(b, i);

      // Make NA if the pixel is invalid
      for(int n = 0; n < 8; n++) {

        // bool out_of_range = R_IsNA(neigh_b[n]); // this did not work
        bool out_of_range = neigh_b[n] == out_value; // it's NA for integers in cpp

        // First check if cell is out of range, because if true, we cant evaluate
        // this cell in the {burning, burnable, burned} vectors
        if(out_of_range) {
          neigh_b[n] = NA_INTEGER;
          positions[n] = NA_INTEGER;
        }
        // If it's in range, evaluate whether it's currently burnable
        else {
          bool burnable_cell = (burning[neigh_b[n]] == 0) &
            (burnable[neigh_b[n]] == 1) &
            (burned[neigh_b[n]] == 0);

          if(!burnable_cell) {
            neigh_b[n] = NA_INTEGER;
            positions[n] = NA_INTEGER;
          }
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
        burning[neigh_b] = spread_around_cpp(
          data_burning,
          data_neighbours,
          coef,
          positions, // columns of valid neighbours
          wind_column,
          elev_column,
          distances,
          upper_limit
        );


      }
    } // end loop across burning cells

    burning[burning == 2] = 0; // turn off previously burned cells
    burned[burning == 1] = 1;  // set as burned the currently burning.

    // Rcout << "burning:\n" << burning << "\n";

    burning_size = sum(burning);
  } // end while

  return burned;
}
