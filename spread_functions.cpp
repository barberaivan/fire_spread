#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <vector>

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
const double elevation_mean = 1163.3;
const double elevation_sd = 399.5;
const int moves[8][2] = {
  {-1, -1},
  {-1,  0},
  {-1,  1},
  { 0, -1},
  { 0,  1},
  { 1, -1},
  { 1,  0},
  { 1,  1}
};

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

const double angles[8] = {
  M_PI * 3 / 4, M_PI, M_PI * 5 / 4,
      M_PI / 2,       M_PI * 3 / 2,
      M_PI / 4,    0, M_PI * 7 / 4
};

// --------------------------------------------------------------------------

//' @title spread_onepix_prob_cpp
//' @description Calculates the probability of a cell spreading fire to another.
//' @return double [0, 1] indicating the probability.
//'
//' @param arma::rowvec data_burning: environmental data from burning cell.
//' @param arma::rowvec data_neighbour: environmental data from target neighbour.
//' @param arma::rowvec coef: parameters in logistic regression to compute the
//'   spread probability as a function of covariates. It has one more elements
//'   than the columns of data_ because it includes the intercept.
//' @param int position: relative position of the neighbour in relation to the
//' burning cell. The eight neighbours are labelled from 0 to 7 beginning from
//' the upper-left one (by row):
//'   0 1 2
//'   3   4
//'   5 6 7.
//'   This is necessary to compute the elevation and wind effects, as they
//'   depend on the angle and distance between burning and target pixels.
//' @param int wind_column: column in the data (landscape) with wind value.
//' @param int elev_column: column in the data (landscape) with elevation value.
//'   Wind and elevation columns must be the last 2.
//' @param double distances: distance (m) between burning and target cell. Used
//' to compute the elevation effect.
//'   This vector depends on the neighbourhood design and on the pixel scale.
//'   If unchanged, it's always the same.
//' @param double upper_limit: upper limit for spread probability (setting to
//'   1 makes absurdly large fires).

// [[Rcpp::export]]
double spread_onepix_prob_cpp(arma::rowvec data_burning,
                           arma::rowvec data_neighbour,
                           arma::rowvec coef,
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
  double wind_term = cos(angles[position] - data_burning(wind_column));

  // elevation term (standardize predictor)
  double elev_term = (data_neighbour(elev_column) - elevation_mean) / elevation_sd;

  // compute probability
  double linpred = coef[0]; // linear predictor, initialize as intercept.
  // (the design matrix lacks the intercept column.)
  // All the effects besides elevation, slope and slope:
  for(int k = 0; k < (data_neighbour.size() - 2); k++) {
    linpred += coef[k+1] * data_neighbour(k);
  }

  // (coef is lagged ahead wrt the data_neighbours because it has the intercept,
  // which is not present in the data.)
  // Add elevation, slope and wind effects
  linpred += wind_term * coef[wind_column + 1] +
    elev_term * coef[elev_column + 1] +
    slope_term * coef[elev_column + 2]; // slope coef is after the elevation one

  double prob = upper_limit / (1 + exp(-linpred));

  return prob;
}



// The same but evaluating the probability (here for backwards compatibility)

// [[Rcpp::export]]
int spread_onepix_cpp(arma::rowvec data_burning,
                      arma::rowvec data_neighbour,
                      arma::rowvec coef,
                      int position,
                      int wind_column,
                      int elev_column,
                      double distance,
                      double upper_limit = 1.0) {

  double prob = spread_onepix_prob_cpp(
    data_burning,
    data_neighbour,
    coef,
    position,
    wind_column,
    elev_column,
    distance,
    upper_limit
  );

  return (int)R::rbinom(1.0, prob);
}

// -----------------------------------------------------------------------

typedef struct _s_burned_res {
  IntegerMatrix burned_bin;
  IntegerMatrix burned_ids;
  int end;
} burned_res;

// function to simulate a fire spread given the landscape,
// model coefficients and ignition points.

burned_res simulate_fire_internal(
    arma::cube landscape,
    IntegerMatrix ignition_cells,
    IntegerMatrix burnable,
    arma::rowvec coef,
    int wind_layer,
    int elev_layer,
    arma::rowvec distances,
    double upper_limit = 1.0,
    double (*prob_fn)(double, double) = R::rbinom) {

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
      arma::rowvec data_burning = landscape.tube(burned_ids(0, b), burned_ids(1, b));

      int neighbours[2][8];
      // get neighbours (adjacent computation here)
      for(int i = 0; i < 8; i++) {
        neighbours[0][i] = burned_ids(0, b) + moves[i][0];
        neighbours[1][i] = burned_ids(1, b) + moves[i][1];
      }


      // Loop over neighbours of the focal burning cell

      for(int n = 0; n < 8; n++) {

        // Is the cell in range?
        bool out_of_range = (
          (neighbours[0][n] < 0) | (neighbours[0][n] >= n_row) | // check rows
          (neighbours[1][n] < 0) | (neighbours[1][n] >= n_col)   // check cols
        );
        if(out_of_range) continue;

        // Is the cell burnable?
        bool burnable_cell = (burned_bin(neighbours[0][n], neighbours[1][n]) == 0) &
                             (burnable(neighbours[0][n], neighbours[1][n]) == 1);

        if(!burnable_cell) continue;

        // obtain data from the neighbour
        arma::rowvec data_neighbour = landscape.tube(neighbours[0][n], neighbours[1][n]);

        // simulate fire
        double prob = spread_onepix_prob_cpp(
          data_burning,
          data_neighbour,
          coef,
          n,           // pixel position identifier (for wind and slope effects)
          wind_layer,
          elev_layer,
          distances[n],
          upper_limit
        );

        int burn = int(prob_fn(1.0, prob));

        if(burn == 0) continue;

        // If burned,
        // store id of recently burned cell and
        // set 1 in burned_bin
        // (but advance end_forward first)
        end_forward += 1;
        burned_ids(0, end_forward) = neighbours[0][n];
        burned_ids(1, end_forward) = neighbours[1][n];
        burned_bin(neighbours[0][n], neighbours[1][n]) = 1;

      } // end loop over neighbours of burning cell b

    } // end loop over burning cells from this cycle

    // update start and end
    start = end + 1;
    end = end_forward;
    burning_size = end - start + 1;

  } // end while

  return {burned_bin, burned_ids, end};
}

//' @title simulate_fire_cpp
//' @description function to simulate a fire spread given the landscape,
//'   model coefficients and ignition points.
//' @return IntegerVector(n_row, n_col): burned layer, coded as
//'   0 = not burned,
//'   1 = burned.
//'   This vector is converted from the burned layer (which is a matrix) by row,
//'   as terra does.

//' @param arma::cube landscape: environmental data from the whole landscape.
//'   See description in spread_around. Every element in the list is a layer in
//'   matrix representation, with rows and columns displayed as they are in the
//'   landscape (not transposed!).
//' @param IntegerMatrix ignition_cells(2, burning_cells): row and column id for
//'   the cell(s) where the fire begun. First row has the row_id, second row has
//'   the col_id.
//' @param IntegerMatrix burnable: matrix indicating if each pixel is burnable (1)
//'   or not (0).
//' @param arma::rowvec coef: parameters in logistic regression to compute the
//'   spread probability as a function of covariates.
//' @param int wind_layer: layer in the data (landscape) with wind matrix.
//' @param int elev_layer: layer in the data (landscape) with elevation matrix.
//'   Wind and elevation layers must be the last 2.
//' @param arma::rowvec distances: distances (m) between burning and target cells,
//'   in the same order as positions. Used to compute the elevation effect.
//'   This vector depends on the neighbourhood design and on the pixel scale.
//'   If unchanged, it's always the same.
//' @param double upper_limit: upper limit for spread probability (setting to
//'   1 makes absurdly large fires).

// [[Rcpp::export]]
IntegerMatrix simulate_fire_cpp(
    arma::cube landscape,
    IntegerMatrix ignition_cells,
    IntegerMatrix burnable,
    arma::rowvec coef,
    int wind_layer,
    int elev_layer,
    arma::rowvec distances,
    double upper_limit) {

  return simulate_fire_internal(
    landscape,
    ignition_cells,
    burnable,
    coef,
    wind_layer,
    elev_layer,
    distances,
    upper_limit
  ).burned_bin;
}



// ...........................................................................
// The same function but deterministic, to test if the discrepancy between R and
// cpp is caused by seed problems

// [[Rcpp::export]]
IntegerMatrix simulate_fire_deterministic_cpp(
    arma::cube landscape,
    IntegerMatrix ignition_cells,
    IntegerMatrix burnable,
    arma::rowvec coef,
    int wind_layer,
    int elev_layer,
    arma::rowvec distances,
    double upper_limit) {
  return simulate_fire_internal(
    landscape,
    ignition_cells,
    burnable,
    coef,
    wind_layer,
    elev_layer,
    distances,
    upper_limit,
    [](double _, double x) { return (double)(x >= 0.5); }
  ).burned_bin;
}

// -------------------------------------------------------------------------

typedef struct _s_burned_compare {
  IntegerMatrix burned_layer;
  IntegerMatrix burned_ids;
  NumericVector counts_veg;
} burned_compare;

// same as simulate_fire, but returning many things to compute discrepancy
// or similarity metrics:
// burned_bin layer,
// burned_ids,
// size by veg_type

burned_compare simulate_fire_compare_cpp(
    arma::cube landscape,
    IntegerMatrix ignition_cells,
    IntegerMatrix burnable,
    arma::rowvec coef,
    int wind_layer,
    int elev_layer,
    arma::rowvec distances,
    double upper_limit) {


  burned_res burned_bin_ids = simulate_fire_internal(
    landscape,
    ignition_cells,
    burnable,
    coef,
    wind_layer,
    elev_layer,
    distances,
    upper_limit
  );

  IntegerMatrix burned_bin = burned_bin_ids.burned_bin;
  IntegerMatrix burned_ids = burned_bin_ids.burned_ids;
  int end = burned_bin_ids.end;


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

  return {burned_bin, burned_ids(_, seq(0, end)), counts_veg};
}

// [[Rcpp::export]]
List simulate_fire_compare(
    arma::cube landscape,
    IntegerMatrix ignition_cells,
    IntegerMatrix burnable,
    arma::rowvec coef,
    int wind_layer,
    int elev_layer,
    arma::rowvec distances,
    double upper_limit) {


  burned_compare burned_com = simulate_fire_compare_cpp(
    landscape,
    ignition_cells,
    burnable,
    coef,
    wind_layer,
    elev_layer,
    distances,
    upper_limit
  );

  // List to return:
  List L = List::create(Named("burned_layer") = burned_com.burned_layer,
                        Named("burned_ids") = burned_com.burned_ids,
                        Named("counts_veg") = burned_com.counts_veg);

  return L;
}