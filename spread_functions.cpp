#include <RcppArmadillo.h>
// armadillo used for the matrix representation of the rasters, to use the
// cube data type.

// Useful armadillo links
//   https://dcgerard.github.io/advancedr/08_cpp_armadillo.html
//   https://arma.sourceforge.net/docs.html#subcube

using namespace Rcpp;

/*
 * Functions to spread fire.
 *
 * The fire spread model is a cellular automata that spreads towards the 8-
 * neighbouring pixels.
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
 * neighbour. They are ordered like this:
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


// --------------------------------------------------------------------------

//' @title spread_onepix_cpp
//' @description Spreads fire from a burning cell to just one neighbour. It's
//'   just a not-vectorized version of spread_around_cpp.
//' @return int burn {0, 1} indicating whether the pixel burned or not.
//'
//' @param arma::vec data_burning: environmental data from burning cell.
//' @param arma::vec data_neighbours: environmental data from target neighbours.
//' @param NumericVector coef: parameters in logistic regression to compute the
//'   spread probability as a function of covariates. It has one more elements
//'   than the columns of data_ because it includes the intercept.
//' @param IntegerVector position: relative position of the neighbour in
//'   relation to the burning cell. The eight neighbours are labelled from 0 to
//'   7 beginning from the upper-left one (by row):
//'   0 1 2
//'   3   4
//'   5 6 7.
//'   This is necessary to compute the elevation and wind effects, as they
//'   depend on the angle and distance between burning and target pixels.
//' @param int wind_column: column in the data (landscape) with wind value.
//' @param int elev_column: column in the data (landscape) with elevation value.
//'   Wind and elevation columns must be the last 2.
//' @param NumericVector distances: distances (m) between burning and target cells,
//'   in the same order as positions. Used to compute the elevation effect.
//'   This vector depends on the neighbourhood design and on the pixel scale.
//'   If unchanged, it's always the same.
//' @param double upper_limit: upper limit for spread probability (setting to
//'   1 makes absurdly large fires).

//' For further comments see documentation for spread_around_cpp.

// [[Rcpp::export]]
int spread_onepix_cpp(arma::vec data_burning,
                      arma::vec data_neighbour,
                      NumericVector coef,
                      int position,
                      int wind_column,
                      int elev_column,
                      double distance,
                      double upper_limit = 1.0) {

  // compute elevation, slope and wind effects

  // slope (from elevation)
  double slope_term = sin(atan(
    (data_neighbour(elev_column) - data_burning(elev_column)) / distance
    ));

  // wind term
  double wind_term = cos(angles(position) - data_burning(wind_column));

  // elevation term (standardize predictor)
  double elev_term = (data_neighbour(elev_column) - elevation_mean) / elevation_sd;

  // compute probability
  double linpred = coef(0); // linear predictor, initialize as intercept.
  // (the design matrix lacks the intercept column.)
  // All the effects besides elevation, slope and slope:
  for(int k = 0; k < (data_neighbour.size() - 2); k++) {
    linpred += coef(k+1) * data_neighbour(k);
  }

  // (coef is lagged ahead wrt the data_neighbours because it has the intercept,
  // which is not present in the data.)
  // Add elevation, slope and wind effects
  linpred += wind_term * coef[wind_column + 1] +
    elev_term * coef[elev_column + 1] +
    slope_term * coef[elev_column + 2]; // slope coef is after the elevation one

  double prob = upper_limit / (1 + exp(-linpred));

  // Simulate the spread
  int burn = R::rbinom(1.0, prob);

  return burn;

}


// The same but returning probability (to test)

// [[Rcpp::export]]
double spread_onepix_prob_cpp(arma::vec data_burning,
                           arma::vec data_neighbour,
                           NumericVector coef,
                           int position,
                           int wind_column,
                           int elev_column,
                           double distance,
                           double upper_limit = 1.0) {

  // compute elevation, slope and wind effects

  // slope (from elevation)
  double slope_term = sin(atan(
    (data_neighbour(elev_column) - data_burning(elev_column)) / distance
  ));

  // wind term
  double wind_term = cos(angles(position) - data_burning(wind_column));

  // elevation term (standardize predictor)
  double elev_term = (data_neighbour(elev_column) - elevation_mean) / elevation_sd;

  // compute probability
  double linpred = coef(0); // linear predictor, initialize as intercept.
  // (the design matrix lacks the intercept column.)
  // All the effects besides elevation, slope and slope:
  for(int k = 0; k < (data_neighbour.size() - 2); k++) {
    linpred += coef(k+1) * data_neighbour(k);
  }

  // (coef is lagged ahead wrt the data_neighbours because it has the intercept,
  // which is not present in the data.)
  // Add elevation, slope and wind effects
  linpred += wind_term * coef[wind_column + 1] +
    elev_term * coef[elev_column + 1] +
    slope_term * coef[elev_column + 2]; // slope coef is after the elevation one

  double prob = upper_limit / (1 + exp(-linpred));

  // Simulate the spread
  int burn = R::rbinom(1.0, prob);

  return prob;

}

// -----------------------------------------------------------------------

// function to simulate a fire spread given the landscape,
// model coefficients and ignition points.


// NOT EDITED, it's not trivial to work with arrays in Rcpp, so I leave it for
// another time.


//' @title simulate_fire_cpp
//' @description function to simulate a fire spread given the landscape,
//'   model coefficients and ignition points.
//' @return IntegerVector(n_row, n_col): burned layer, coded as
//'   0 = not burned,
//'   1 = burned.
//'   This vector is converted from the burned layer (which is a matrix) by row,
//'   as terra does.

//' @param List landscape: environmental data from the whole landscape.
//'   See description in spread_around. Every element in the list is a layer in
//'   matrix representation, with rows and columns displayed as they are in the
//'   landscape (not transposed!).
//' @param IntegerMatrix ignition_cells(2, burning_cells): row and column id for
//'   the cell(s) where the fire begun. First row has the row_id, second row has
//'   the col_id.
//' @param IntegerMatrix burnable: matrix indicating if each pixel is burnable (1)
//'   or not (0).
//' @param NumericVector coef: parameters in logistic regression to compute the
//'   spread probability as a function of covariates.
//' @param IntegerVector positions: position of each neighbour in
//'   relation to the burning cell. The eight neighbours are labelled from 0 to
//'   7 beginning from the upper-left corner (by row):
//'   0 1 2
//'   3   4
//'   5 6 7.
//'   This is necessary to compute the elevation and wind effects, as they
//'   depend on the angle and distance between burning and target pixels.
//'   (Sometimes not all neighbours are burnable, so this vector will not
//'   always be complete.)
//' @param IntegerMatrix moves: see adjacent_cpp0 documentation.
//' @param int wind_layer: layer in the data (landscape) with wind matrix.
//' @param int elev_layer: layer in the data (landscape) with elevation matrix.
//'   Wind and elevation layers must be the last 2.
//' @param NumericVector distances: distances (m) between burning and target cells,
//'   in the same order as positions. Used to compute the elevation effect.
//'   This vector depends on the neighbourhood design and on the pixel scale.
//'   If unchanged, it's always the same.
//' @param double upper_limit: upper limit for spread probability (setting to
//'   1 makes absurdly large fires).


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
IntegerMatrix simulate_fire_cpp(
    arma::cube landscape,
    IntegerMatrix ignition_cells,
    IntegerMatrix burnable,
    NumericVector coef,
    int wind_layer,
    int elev_layer,
    NumericVector distances,
    double upper_limit) {

  int n_row = burnable.nrow();
  int n_col = burnable.ncol();
  int n_cell = n_row * n_col;

  // burned_ids [row-col, cell] will be filled with the row_col ids (rows) of the
  // burning pixels (columns). start and end integers will define the positions
  // limits corresponding to the burning cells in every burn cycle.
  IntegerMatrix burned_ids(2, n_cell); // check it's filled with 0 // -2147483648 is NA_INTEGER

  int start = 0;
  // end is the last non-empty position in the burned_ids matrix.
  int end = ignition_cells.ncol() - 1;
  // Example:
  // burned_ids = {231, 455, 342, 243, NA, NA, NA, NA};
  //               start          end.
  // if only one cell is burning, start = end.

  // initialize burned_ids and burning_size with ignition_cells
  for(int c = 0; c <= end; c++) {
    for(int r = 0; r < 2; r++) {
      burned_ids(r, c) = ignition_cells(r, c);
    }
  }

  // initialize burning_size
  int burning_size = ignition_cells.ncol(); // == end + 1 - start

  // The burned_bin matrix will indicate whether each pixel is burned or burning
  // (1) or not (0). It's necessary to have this now because it will be used
  // to define burnable neighbours.
  IntegerMatrix burned_bin(n_row, n_col);

  // initialize with ignition_cells
  for(int i = 0; i <= end; i++) {
    burned_bin(ignition_cells(0, i), ignition_cells(1, i)) = 1;
  }

  while(burning_size > 0) {
    // Loop over all the burning cells to burn their neighbours. Use end_forward
    // to update the last position in burned_ids within this loop, without
    // compromising the loop's integrity.
    int end_forward = end;

    // Loop over burning cells in the cycle

    // b is going to keep the position in burned_ids that have to be evaluated
    // in this burn cycle
    for(int b = start; b <= end; b++) {

      // Get burning_cells' data
      arma::vec data_burning = landscape.tube(burned_ids(0, b), burned_ids(1, b));

      // get neighbours (adjacent computation here)
      IntegerMatrix neighbours(2, 8);
      for(int i = 0; i < 8; i++) neighbours(_, i) = burned_ids(_, b) + moves(_, i);

      // Loop over neighbours of the focal burning cell

      for(int n = 0; n < 8; n++) {

        // Is the cell in range?
        bool out_of_range = (
          (neighbours(0, n) < 0) | (neighbours(0, n) >= n_row) | // check rows
          (neighbours(1, n) < 0) | (neighbours(1, n) >= n_col)   // check cols
        );
        if(out_of_range) continue;

        // Is the cell burnable?
        bool burnable_cell = (burned_bin(neighbours(0, n), neighbours(1, n)) == 0) &
                             (burnable(neighbours(0, n), neighbours(1, n)) == 1);

        if(!burnable_cell) continue;

        // obtain data from the neighbour
        arma::vec data_neighbour = landscape.tube(neighbours(0, n), neighbours(1, n));

        // simulate fire
        int burn;
        burn = spread_onepix_cpp(
          data_burning,
          data_neighbour,
          coef,
          n,           // pixel position identifier (for wind and slope effects)
          wind_layer,
          elev_layer,
          distances(n),
          upper_limit
        );

        if(burn == 0) continue;

        // If burned,
        // store id of recently burned cell and
        // set 1 in burned_bin
        // (but advance end_forward first)
        end_forward += 1;
        burned_ids(0, end_forward) = neighbours(0, n);
        burned_ids(1, end_forward) = neighbours(1, n);
        burned_bin(neighbours(0, n), neighbours(1, n)) = 1;

      } // end loop over neighbours of burning cell b

    } // end loop over burning cells from this cycle

    // update start and end
    start = end + 1;
    end = end_forward;
    burning_size = end - start + 1;

  } // end while

  return burned_bin;
}



// ...........................................................................
// The same function but deterministic, to test if the discrepancy between R and
// cpp is caused by seed problems

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
IntegerMatrix simulate_fire_deterministic_cpp(
    arma::cube landscape,
    IntegerMatrix ignition_cells,
    IntegerMatrix burnable,
    NumericVector coef,
    int wind_layer,
    int elev_layer,
    NumericVector distances,
    double upper_limit) {

  int n_row = burnable.nrow();
  int n_col = burnable.ncol();
  int n_cell = n_row * n_col;

  // burned_ids [row-col, cell] will be filled with the row_col ids (rows) of the
  // burning pixels (columns). start and end integers will define the positions
  // limits corresponding to the burning cells in every burn cycle.
  IntegerMatrix burned_ids(2, n_cell); // check it's filled with 0 // -2147483648 is NA_INTEGER

  int start = 0;
  // end is the last non-empty position in the burned_ids matrix.
  int end = ignition_cells.ncol() - 1;
  // Example:
  // burned_ids = {231, 455, 342, 243, NA, NA, NA, NA};
  //               start          end.
  // if only one cell is burning, start = end.

  // initialize burned_ids and burning_size with ignition_cells
  for(int c = 0; c <= end; c++) {
    for(int r = 0; r < 2; r++) {
      burned_ids(r, c) = ignition_cells(r, c);
    }
  }

  // initialize burning_size
  int burning_size = ignition_cells.ncol(); // == end + 1 - start

  // The burned_bin matrix will indicate whether each pixel is burned or burning
  // (1) or not (0). It's necessary to have this now because it will be used
  // to define burnable neighbours.
  IntegerMatrix burned_bin(n_row, n_col);

  // initialize with ignition_cells
  for(int i = 0; i <= end; i++) {
    burned_bin(ignition_cells(0, i), ignition_cells(1, i)) = 1;
  }

  while(burning_size > 0) {
    // Loop over all the burning cells to burn their neighbours. Use end_forward
    // to update the last position in burned_ids within this loop, without
    // compromising the loop's integrity.
    int end_forward = end;

    // Loop over burning cells in the cycle

    // b is going to keep the position in burned_ids that have to be evaluated
    // in this burn cycle
    for(int b = start; b <= end; b++) {

      // Get burning_cells' data
      arma::vec data_burning = landscape.tube(burned_ids(0, b), burned_ids(1, b));

      // get neighbours (adjacent computation here)
      IntegerMatrix neighbours(2, 8);
      for(int i = 0; i < 8; i++) neighbours(_, i) = burned_ids(_, b) + moves(_, i);

      // Loop over neighbours of the focal burning cell

      for(int n = 0; n < 8; n++) {

        // Is the cell in range?
        bool out_of_range = (
          (neighbours(0, n) < 0) | (neighbours(0, n) >= n_row) | // check rows
          (neighbours(1, n) < 0) | (neighbours(1, n) >= n_col)   // check cols
        );
        if(out_of_range) continue;

        // Is the cell burnable?
        bool burnable_cell = (burned_bin(neighbours(0, n), neighbours(1, n)) == 0) &
                             (burnable(neighbours(0, n), neighbours(1, n)) == 1);

        if(!burnable_cell) continue;

        // obtain data from the neighbour
        arma::vec data_neighbour = landscape.tube(neighbours(0, n), neighbours(1, n));

        // simulate fire
        int burn;
        double pburn = spread_onepix_prob_cpp(
          data_burning,
          data_neighbour,
          coef,
          n,           // pixel position identifier (for wind and slope effects)
          wind_layer,
          elev_layer,
          distances(n),
          upper_limit
        );

        IntegerVector fc = neighbours(_, n);

        burn = R::rbinom(1, pburn);

        //// make deterministic here
        if(pburn < 0.5000000000) continue;

        // If burned,
        // store id of recently burned cell and
        // set 1 in burned_bin
        // (but advance end_forward first)
        end_forward += 1;
        burned_ids(0, end_forward) = neighbours(0, n);
        burned_ids(1, end_forward) = neighbours(1, n);
        burned_bin(neighbours(0, n), neighbours(1, n)) = 1;
      } // end loop over neighbours of burning cell b

    } // end loop over burning cells from this cycle

    // update start and end
    start = end + 1;
    end = end_forward;
    burning_size = end - start + 1;

  } // end while

  return burned_bin;
}

// -------------------------------------------------------------------------

// same as simulate_fire, but returning many things to compute discrepancy
// or similarity metrics:
// burned_bin layer,
// burned_ids,
// size by veg_type

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List simulate_fire_compare(
    arma::cube landscape,
    IntegerMatrix ignition_cells,
    IntegerMatrix burnable,
    NumericVector coef,
    int wind_layer,
    int elev_layer,
    NumericVector distances,
    double upper_limit) {

  int n_row = burnable.nrow();
  int n_col = burnable.ncol();
  int n_cell = n_row * n_col;

  // burned_ids [row-col, cell] will be filled with the row_col ids (rows) of the
  // burning pixels (columns). start and end integers will define the positions
  // limits corresponding to the burning cells in every burn cycle.
  IntegerMatrix burned_ids(2, n_cell); // check it's filled with 0 // -2147483648 is NA_INTEGER

  int start = 0;
  // end is the last non-empty position in the burned_ids matrix.
  int end = ignition_cells.ncol() - 1;
  // Example:
  // burned_ids = {231, 455, 342, 243, NA, NA, NA, NA};
  //               start          end.
  // if only one cell is burning, start = end.

  // initialize burned_ids and burning_size with ignition_cells
  for(int c = 0; c <= end; c++) {
    for(int r = 0; r < 2; r++) {
      burned_ids(r, c) = ignition_cells(r, c);
    }
  }

  // initialize burning_size
  int burning_size = ignition_cells.ncol(); // == end + 1 - start

  // The burned_bin matrix will indicate whether each pixel is burned or burning
  // (1) or not (0). It's necessary to have this now because it will be used
  // to define burnable neighbours.
  IntegerMatrix burned_bin(n_row, n_col);

  // initialize with ignition_cells
  for(int i = 0; i <= end; i++) {
    burned_bin(ignition_cells(0, i), ignition_cells(1, i)) = 1;
  }

  while(burning_size > 0) {
    // Loop over all the burning cells to burn their neighbours. Use end_forward
    // to update the last position in burned_ids within this loop, without
    // compromising the loop's integrity.
    int end_forward = end;

    // Loop over burning cells in the cycle

    // b is going to keep the position in burned_ids that have to be evaluated
    // in this burn cycle
    for(int b = start; b <= end; b++) {

      // Get burning_cells' data
      arma::vec data_burning = landscape.tube(burned_ids(0, b), burned_ids(1, b));

      // get neighbours (adjacent computation here)
      IntegerMatrix neighbours(2, 8);
      for(int i = 0; i < 8; i++) neighbours(_, i) = burned_ids(_, b) + moves(_, i);

      // Loop over neighbours of the focal burning cell

      for(int n = 0; n < 8; n++) {

        // Is the cell in range?
        bool out_of_range = (
          (neighbours(0, n) < 0) | (neighbours(0, n) >= n_row) | // check rows
            (neighbours(1, n) < 0) | (neighbours(1, n) >= n_col)   // check cols
        );
        if(out_of_range) continue;

        // Is the cell burnable?
        bool burnable_cell = (burned_bin(neighbours(0, n), neighbours(1, n)) == 0) &
          (burnable(neighbours(0, n), neighbours(1, n)) == 1);

        if(!burnable_cell) continue;

        // obtain data from the neighbour
        arma::vec data_neighbour = landscape.tube(neighbours(0, n), neighbours(1, n));

        // simulate fire
        int burn;
        burn = spread_onepix_cpp(
          data_burning,
          data_neighbour,
          coef,
          n,           // pixel position identifier (for wind and slope effects)
          wind_layer,
          elev_layer,
          distances(n),
          upper_limit
        );

        if(burn == 0) continue;

        // If burned,
        // store id of recently burned cell and
        // set 1 in burned_bin
        // (but advance end_forward first)
        end_forward += 1;
        burned_ids(0, end_forward) = neighbours(0, n);
        burned_ids(1, end_forward) = neighbours(1, n);
        burned_bin(neighbours(0, n), neighbours(1, n)) = 1;

      } // end loop over neighbours of burning cell b

    } // end loop over burning cells from this cycle

    // update start and end
    start = end + 1;
    end = end_forward;
    burning_size = end - start + 1;

  } // end while

  // Compute burned area by vegetation type
  NumericVector counts_veg(4);
  for(int i = 0; i <= end; i++) {
    for(int v = 1; v < 4; v++) { // starts from 1 because 0 is for shrubland
      counts_veg(v) += landscape(burned_ids(0, i), burned_ids(1, i), v-1);
    }
    // landscape must have the vegetation layers in the first veg_types - 1 layers
  }
  // fill shrubland
  counts_veg(0) = (end + 1) - sum(counts_veg[seq(1, 3)]);

  // List to return:
  List L = List::create(Named("burned_layer") = burned_bin,
                        Named("burned_ids") = burned_ids(_, seq(0, end)),
                        Named("counts_veg") = counts_veg);

  return L;
}


// simulate_fire but not computing prob -------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
IntegerMatrix simulate_fire_notprob(
    arma::cube landscape,
    IntegerMatrix ignition_cells,
    IntegerMatrix burnable,
    NumericVector coef,
    int wind_layer,
    int elev_layer,
    NumericVector distances,
    double upper_limit) {
  
  int n_row = burnable.nrow();
  int n_col = burnable.ncol();
  int n_cell = n_row * n_col;
  
  // burned_ids [row-col, cell] will be filled with the row_col ids (rows) of the
  // burning pixels (columns). start and end integers will define the positions
  // limits corresponding to the burning cells in every burn cycle.
  IntegerMatrix burned_ids(2, n_cell); // check it's filled with 0 // -2147483648 is NA_INTEGER
  
  int start = 0;
  // end is the last non-empty position in the burned_ids matrix.
  int end = ignition_cells.ncol() - 1;
  // Example:
  // burned_ids = {231, 455, 342, 243, NA, NA, NA, NA};
  //               start          end.
  // if only one cell is burning, start = end.
  
  // initialize burned_ids and burning_size with ignition_cells
  for(int c = 0; c <= end; c++) {
    for(int r = 0; r < 2; r++) {
      burned_ids(r, c) = ignition_cells(r, c);
    }
  }
  
  // initialize burning_size
  int burning_size = ignition_cells.ncol(); // == end + 1 - start
  
  // The burned_bin matrix will indicate whether each pixel is burned or burning
  // (1) or not (0). It's necessary to have this now because it will be used
  // to define burnable neighbours.
  IntegerMatrix burned_bin(n_row, n_col);
  
  // initialize with ignition_cells
  for(int i = 0; i <= end; i++) {
    burned_bin(ignition_cells(0, i), ignition_cells(1, i)) = 1;
  }
  
  while(burning_size > 0) {
    // Loop over all the burning cells to burn their neighbours. Use end_forward
    // to update the last position in burned_ids within this loop, without
    // compromising the loop's integrity.
    int end_forward = end;
    
    // Loop over burning cells in the cycle
    
    // b is going to keep the position in burned_ids that have to be evaluated
    // in this burn cycle
    for(int b = start; b <= end; b++) {
      
      // Get burning_cells' data
      arma::vec data_burning = landscape.tube(burned_ids(0, b), burned_ids(1, b));
      
      // get neighbours (adjacent computation here)
      IntegerMatrix neighbours(2, 8);
      for(int i = 0; i < 8; i++) neighbours(_, i) = burned_ids(_, b) + moves(_, i);
      
      // Loop over neighbours of the focal burning cell
      
      for(int n = 0; n < 8; n++) {
        
        // Is the cell in range?
        bool out_of_range = (
          (neighbours(0, n) < 0) | (neighbours(0, n) >= n_row) | // check rows
            (neighbours(1, n) < 0) | (neighbours(1, n) >= n_col)   // check cols
        );
        if(out_of_range) continue;
        
        // Is the cell burnable?
        bool burnable_cell = (burned_bin(neighbours(0, n), neighbours(1, n)) == 0) &
          (burnable(neighbours(0, n), neighbours(1, n)) == 1);
        
        if(!burnable_cell) continue;
        
        // obtain data from the neighbour
        arma::vec data_neighbour = landscape.tube(neighbours(0, n), neighbours(1, n));
        
        // simulate fire
        int burn;

        double pburn = landscape(1, 1, 1); // subset to mimic the extraction of the saved prob
        pburn = 1.0;
        // double pburn = spread_onepix_prob_cpp(
        //   data_burning,
        //   data_neighbour,
        //   coef,
        //   n,           // pixel position identifier (for wind and slope effects)
        //   wind_layer,
        //   elev_layer,
        //   distances(n),
        //   upper_limit
        // );
        
        burn = R::rbinom(1, pburn);
        
        //// make deterministic here
        if(pburn < 0.5000000000) continue;
        
        // If burned,
        // store id of recently burned cell and
        // set 1 in burned_bin
        // (but advance end_forward first)
        end_forward += 1;
        burned_ids(0, end_forward) = neighbours(0, n);
        burned_ids(1, end_forward) = neighbours(1, n);
        burned_bin(neighbours(0, n), neighbours(1, n)) = 1;
        
      } // end loop over neighbours of burning cell b
      
    } // end loop over burning cells from this cycle
    
    // update start and end
    start = end + 1;
    end = end_forward;
    burning_size = end - start + 1;
    
  } // end while
  
  return burned_bin;
}

// -----------------------------------------------------------------------


//' @title compare_fires_try  
//' @description Function compare two fires using many similarity
//'   indexes, to try them as proxies for the likelihood, when simulated fires 
//'   are compared to the observed one. "_try" because it computes many 
//'   similarity indexes; after selecting one we will have a function to compute
//'   only the best one.
//' @return NumericVector(n_metrics): vector with the comparison indexes.

//' @param List fire1, List fire2: data from the fires to compare. This has the
//'   same elements as the result from simulate_fire_compare: 
//'     IntegerMatrix burned_layer: binary matrix storing (1) in burned pixels;
//'     IntegerMatrix burned_ids: id in [row, col] (0-indexing) of the burned 
//'       cells. First row holds the rows, second row holds the columns, and 
//'       each column is a burned pixel;
//'     IntegerVector counts_veg: count of burned pixels by vegetation type.
//'   The burned_layer is used to compute the spatial overlap index; the 
//'   burned_ids is used to evaluate the common burned pixels, looping only in 
//'   the burned_ids from the smaller fire; counts_veg is used to compute the 
//'   difference in number of pixels burned by vegetation type.
//' @param double lscale: length-scale parameter for the gaussian kernel 
//'   (the sd in a Normal distribution). Used to turn a dissimilarity into 
//'   a similarity.

//' Details: the discrepancy index computed in Morales et al. 2015 paper is not
//' considered here because it has a wayward behaviour. Different denominators
//' are used, which work better.

// [[Rcpp::export]]
NumericVector compare_fires_try(List fire1, List fire2,
                                double lscale = 0.2) {
  
  // Extract list elemnts ------------------------------------------------
  
  NumericMatrix burned1 = fire1["burned_layer"];
  NumericMatrix burned2 = fire2["burned_layer"];
  
  IntegerMatrix burned_ids1 = fire1["burned_ids"];
  IntegerMatrix burned_ids2 = fire2["burned_ids"];
  
  double size1 = burned_ids1.ncol();
  double size2 = burned_ids2.ncol();
  
  NumericVector counts1 = fire1["counts_veg"];
  NumericVector counts2 = fire2["counts_veg"];
  
  // overlap_sp -----------------------------------------------------------
  
  double common = 0.0;
  // compute common pixels only in the smaller fire
  if(size1 <= size2) {
    for(int i = 0; i < size1; i++) {
      common += burned2(burned_ids1(0, i), burned_ids1(1, i));
    }
  } else {
    for(int i = 0; i < size2; i++) {
      common += burned1(burned_ids2(0, i), burned_ids2(1, i));
    }
  }
  
  double overlap_sp = common / (size1 + size2 - common);
  
  // overlap_vd -----------------------------------------------------------
  
  // Get vegetation distribution by fire (normalized burned areas)
  int veg_types = counts1.length();
  
  NumericVector burned_dist_1(veg_types);
  NumericVector burned_dist_2(veg_types);
  
  for(int v = 0; v < veg_types; v++) {
    burned_dist_1[v] = counts1[v] / sum(counts1);
    burned_dist_2[v] = counts2[v] / sum(counts2);
  }
  
  // compute vegetation distribution overlap
  double overlap_vd = 0.0;
  for(int v = 0; v < veg_types; v++) {
    overlap_vd += std::min(burned_dist_1[v], burned_dist_2[v]);
  }
  
  // deltas by veg_type ---------------------------------------------------
  
  // normalized difference using absolute difference. The difference by veg_type
  // is in [0, 1]. So, if we divide delta_norm by veg_num, it will be in [0, 1].
  double delta_norm = 0.0;
  double delta_pow = 0.0;
  
  for(int v = 0; v < veg_types; v++) {
    double sum_area = counts1[v] + counts2[v];
    if(sum_area > 0.0) {
      delta_norm += std::abs((counts1(v) - counts2(v)) / sum_area);
    }
  }
  
  double veg_num = counts1.length();
  
  // Scale to [0, 1]
  double delta_norm_unit = delta_norm / veg_num;
  
  // Transform to similarities
  double overlap_norm = 1.0 - delta_norm_unit;
  double overlap_expquad = exp(-pow(delta_norm_unit, 2.0) / lscale); // 0.2 is the Gaussian SD.
  double overlap_quad = 1 - pow(delta_norm_unit, 2.0);
  
  // ---------------------------------------------------------------------
  
  NumericVector indexes = NumericVector::create(
    // pure indices
    Named("overlap_sp")      = overlap_sp,
    
    Named("overlap_vd")      = overlap_vd,
    Named("overlap_norm")    = overlap_norm,
    Named("overlap_expquad") = overlap_expquad,
    Named("overlap_quad")    = overlap_quad,
    
    // mixture indices
    Named("sp_norm_5050")    = 0.50 * overlap_sp + 0.50 * overlap_norm,
    Named("sp_norm_7525")    = 0.75 * overlap_sp + 0.25 * overlap_norm,
    Named("sp_expquad_5050") = 0.50 * overlap_sp + 0.50 * overlap_expquad,
    Named("sp_expquad_7525") = 0.75 * overlap_sp + 0.25 * overlap_expquad,
    Named("sp_quad_5050")    = 0.50 * overlap_sp + 0.50 * overlap_quad,
    Named("sp_quad_7525")    = 0.75 * overlap_sp + 0.25 * overlap_quad  
  );
  
  return indexes;
}

// Emulate likelihood function --------------------------------------------
// trying many discrepancy functions

//' @title emulate_loglik_try
//' @description Function to emulate the spread model's likelihood, approximated
//'   by many similarity indexes between the observed and 
//'   simulated fires. It's a wrapper around simulate_fire_compare, which returns
//'   the simulated fire with data useful for comparison. In addition, it takes
//'   as arguments data from the observed fire to be compared. 
//'   "_try" because it computes many similarity indexes; after selecting one 
//'   we will have a function to consider only the best one.
//' @return NumericMatrix(n_replicates, n_metrics): Matrix with each comparison 
//'   metric (colummns) by simulated fire (rows)

//' @param array(3D) landscape: Environmental data from the whole landscape.
//'   See description in spread_around. The 3rd dimension contains each layer
//'   (covariates). The landscape is represented in matrix form.
//' @param IntegerMatrix ignition_cells(2, burning_cells): row and column id for
//'   the cell(s) where the fire begun. First row has the row_id, second row has
//'   the col_id.
//' @param IntegerMatrix burnable: matrix indicating if each pixel is burnable (1)
//'   or not (0).
//' @param NumericVector coef: parameters in logistic regression to compute the
//'   spread probability as a function of covariates.
//' @param int wind_layer: layer in the data (landscape) with wind matrix.
//' @param int elev_layer: layer in the data (landscape) with elevation matrix.
//'   Wind and elevation layers must be the last 2.
//' @param NumericVector distances: distances (m) between burning and target cells,
//'   in the same order as positions. Used to compute the elevation effect.
//'   This vector depends on the neighbourhood design and on the pixel scale.
//' @param double upper_limit: upper limit for spread probability (setting to
//'   1 makes absurdly large fires).
//' @param List fire_obs: Data of the observed (reference) fire. This has the 
//'   same elements as the result from simulate_fire_compare: 
//'     IntegerMatrix burned_layer: binary matrix storing (1) in burned pixels;
//'     IntegerMatrix burned_ids: id in [row, col] of the burned cells. First 
//'       row hols the rows, second row holds the columns, and each column is 
//'       a pixel;
//'     IntegerVector counts_veg: count of burned pixels by vegetation type.
//'   The burned_layer is used to compute the spatial overlap index; the 
//'   burned_ids is used to evaluate the common burned pixels only in the smallest
//'   of the fires, evaluating only those burned_ids; counts_veg is used to
//'   compute the difference in number of pixels burned by vegetation type.
//' @param int n_replicates: number of fires to simulate, defaults to 10.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix emulate_loglik_try(
    // arguments for simulate_fire
    arma::cube landscape,
    IntegerMatrix ignition_cells,
    IntegerMatrix burnable,
    NumericVector coef,
    int wind_layer,
    int elev_layer,
    NumericVector distances,
    double upper_limit,
    // arguments to compare simulated and observed fire
    List fire_ref,
    // number of simulated fires
    int n_replicates = 10,
    // number of indices to compare fires
    int n_indices = 11
) {

  NumericMatrix similarity(n_replicates, n_indices);

  for(int i = 0; i < n_replicates; i++) {

    // simulate_fire
    List fire_sim = simulate_fire_compare(
      landscape,
      ignition_cells,
      burnable,
      coef,
      wind_layer,
      elev_layer,
      distances,
      upper_limit
    );

    similarity(i, _) = compare_fires_try(fire_ref, fire_sim);
  }

  return similarity;
}