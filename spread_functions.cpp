#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <vector>

// armadillo is used for the matrix representation of the rasters, to use the
// fcube data type.

// Useful armadillo links
//   https://dcgerard.github.io/advancedr/08_cpp_armadillo.html
//   https://arma.sourceforge.net/docs.html#subfcube

using namespace Rcpp;

/*
 * Functions to spread fire.
 *
 * The fire spread model is a cellular automata that spreads towards the 8-
 * neighbouring pixels.
 * 
 * The landscape is defined by a vegetation (integer) matrix, with classes from
 * 0 to n_veg_types - 1 (99 is non-burnable) and a terrain arma::fcube 
 * (a 3D array), where each layer is a matrix: 
 * {northing, elevation, wind direction}. The fourth terrain
 * variable is the slope index, which is directional, so it's computed during 
 * the simulation. 
 * 
 * In the most complex desired model n_veg_types = 6 (see below). Landscapes 
 * are produced with 6 classes, and for simpler models with aggregated classes
 * the same parameter is assigned to equivalent classes.
 *  
 * The fire-parameters (coef) have no global intercept term, i.e., each 
 * vegetation type has its own intercept. This is aimed to reduce the 
 * correlation in the likelihood and to facilitate the exclusion of vegetation-
 * parameters corresponding to types absent in the landscape. (Estimating an
 * intercept when the correponding vegetation class is absent is troublesome.)
 * 
 * The coef vector holds all the parameters to run the simulation. The first
 * n_veg_types values are the intercepts for each vegetation type; the remaining 
 * 4 are the terrain parameters (northing, elevation, wind direction, slope). 
 * This vecto is separated as 
 * coef_veg = coef[0 : (n_veg_types - 1)]
 * coef_terrain = coef[n_veg_types : coef.length()]
 * in the simulate_fire() function.
 * 
 * The interannual climatic variability is represented through the 
 * summer-average (dec-march) FWI anomaly at the pixel level, and is treated as
 * constant within a fire or landscape. Its value at the ignition point 
 * is used to define the mean of a fire-level random effect, with a linear 
 * function. This is not included as a pixel-level variable because of its low
 * resolution, which made its values almost constant within landscapes, which
 * would generete high correlation between the fwi and vegetation parameters.
 */

// Constants ---------------------------------------------------------------

// Elevation data to standardize distance between pixels
const float elevation_mean = 1163.3;
const float elevation_sd = 399.5;

// Distance between pixels, in m / elevation_sd.
// 30 m is the landscape resolution. 
const float distances[8] = {
  30 * sqrtf(2) / elevation_sd, 30 / elevation_sd, 30 * sqrtf(2) / elevation_sd,
  30 / elevation_sd,                               30 / elevation_sd,
  30 * sqrtf(2) / elevation_sd, 30 / elevation_sd, 30 * sqrtf(2) / elevation_sd
};
// it didn't let me loop to divide by elevation_sd

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

const float angles[8] = {
  M_PI * 3 / 4, M_PI, M_PI * 5 / 4,
  M_PI / 2,           M_PI * 3 / 2,
  M_PI / 4,      0,   M_PI * 7 / 4
};


// Terrain coefficients names
enum terrain_names {
  northing,
  elev, 
  windir,
  slope
};

// Vegetation coefficients names
enum veg_names {
  shrubland,
  subalpine, 
  wet,
  dry_a, // araucaria
  dry_b, // cypress
  steppe
};

// --------------------------------------------------------------------------

//' @title spread_onepix_prob_cpp
//' @description Calculates the probability of a cell spreading fire to another.
//' @return float [0, 1] indicating the probability.
//'
//' @param int vegetation_type: vegetation type of the target cell.
//' @param arma::frowvec terrain_burning: terrain data from burning cell. 
//' @param arma::frowvec terrain_neighbour: terrain data from target neighbour.
//' @param arma::frowvec coef_veg: intercepts in logistic regression associated
//'   to each vegetation type. To aggregate vegetation types (e.g., dry_forest =
//'   dry_forest_a or dry_forest_b), the same values are assigned to the 
//'   aggregated classes.
//' @param arma::frowvec coef_terrain: slopes in logistic regression related
//'   multiplying the terrain predictors. 
//' @param int position: relative position of the neighbour in relation to the
//' burning cell. The eight neighbours are labelled from 0 to 7 beginning from
//' the upper-left one (by row):
//'   0 1 2
//'   3   4
//'   5 6 7.
//'   This is necessary to compute the slope and wind effects, as they
//'   depend on the angle and distance between burning and target pixels.
//' @param float upper_limit: upper limit for spread probability (setting to
//'   1 makes absurdly large fires; 0.5 is preferred).

// [[Rcpp::export]]
float spread_onepix_prob_cpp(
  int vegetation_type,
  arma::frowvec terrain_burning,
  arma::frowvec terrain_neighbour,
  arma::frowvec coef_veg,
  arma::frowvec coef_terrain,
  int position,
  float upper_limit = 1.0
) {
  
  // wind term
  float wind_term = cosf(angles[position] - terrain_burning(windir));
  
  // slope term (from elevation and distance)
  float slope_term = sinf(atanf(
    (terrain_neighbour(elev) - terrain_burning(elev)) / distances[position]
  ));

  // compute linear predictor
  float linpred = coef_veg[vegetation_type]; // initialize at the corresponding 
                                             // intercept.

  linpred += coef_terrain[northing] * terrain_neighbour[northing] +
             coef_terrain[elev] * terrain_neighbour[elev] +
             coef_terrain[windir] * wind_term +
             coef_terrain[slope] * slope_term;
  
  // burn probability
  float prob = upper_limit / (1 + expf(-linpred));

  return prob;
}

// The same but evaluating the probability and simulating the burn from a
// Bernoulli distribution (here for backwards compatibility)

// [[Rcpp::export]]
int spread_onepix_cpp(
  int vegetation_type,
  arma::frowvec terrain_burning,
  arma::frowvec terrain_neighbour,
  arma::frowvec coef_veg,
  arma::frowvec coef_terrain,
  int position,
  float upper_limit = 1.0
) {

  float prob = spread_onepix_prob_cpp(
    vegetation_type,
    terrain_burning,
    terrain_neighbour,
    coef_veg,
    coef_terrain,
    position,
    upper_limit
  );

  return (int)R::rbinom(1.0, prob);
}

// -----------------------------------------------------------------------

//' @title simulate_fire_internal
//' @description function to simulate a fire spread given the landscape,
//'   model coefficients and ignition points. 
//' @return burned_res: struct containing the following objects:
//'   IntegerMatrix burned_bin, a binary matrix indicating the burned pixels
//'   IntegerMatrix burned_ids, a matrix with a column by burned pixel, 
//'     indicating its row (row1) and column (row2) in the landscape,
//'   int end, the number of burned pixels.

//' @param IntegerMatrix vegetation: integer matrix representing the vegetation 
//'   type. 99 is non-burnable, and valid values are {0, ..., n_veg_types - 1}.
//' @param arma::fcube terrain: terrain data, where each matrix slice is a 
//'   predictor: {northing, elevation, wind direction}. Slope is absent
//'   because it's directional, so it's computed during the simulation.
//' @param int n_veg_types: integer indicating the number of vegetation types
//'   considered by the model (not just those present in the focal landscape).
//'   used to read properly the vegetation- and non-vegetation parameters in 
//'   coef.  
//'   
//' @param IntegerMatrix ignition_cells(2, burning_cells): row and column id for
//'   the cell(s) where the fire begun. First row has the row_id, second row has
//'   the col_id.
//' @param arma::frowvec coef: parameters in logistic regression to compute the
//'   spread probability as a function of covariates.
//' @param float upper_limit: upper limit for spread probability (setting to
//'   1 makes absurdly large fires).
//' @param function [unnamed]: function evaluating the burn probability. By 
//'   default, it's R::rbinom(), but can be set to a deterministic behaviour to 
//'   test whether the probability computation matches R's function. (Just for 
//'   testing.)

typedef struct _s_burned_res {
  IntegerMatrix burned_bin;
  IntegerMatrix burned_ids;
  int end;
} burned_res;

burned_res simulate_fire_internal(
    IntegerMatrix vegetation,
    arma::fcube terrain,
    IntegerMatrix ignition_cells,
    arma::frowvec coef,
    int n_veg_types = 6,
    float upper_limit = 1.0,
    double (*prob_fn)(double, double) = R::rbinom
  ) {
  
  // separate vegetation and terrain coefficients
  arma::frowvec coef_veg(n_veg_types);
  arma::frowvec coef_terrain(4);
  coef_veg = coef.subvec(0, n_veg_types - 1);
  coef_terrain = coef.subvec(n_veg_types, coef.size() - 1);
  
  // define landscape dimensions
  int n_row = vegetation.nrow();
  int n_col = vegetation.ncol();
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

      // Get burning_cell's data
      arma::frowvec terrain_burning = terrain.tube(burned_ids(0, b), burned_ids(1, b));

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
        
        // Get vegetation class to know whether the cell is burnable
        int veg_target = vegetation(neighbours[0][n], neighbours[1][n]);
        
        // Is the cell burnable?
        bool burnable_cell = 
          (burned_bin(neighbours[0][n], neighbours[1][n]) == 0) & // not burned
          (veg_target < 99);                                      // burnable
        if(!burnable_cell) continue;

        // obtain data from the neighbour
        arma::frowvec terrain_neighbour = terrain.tube(neighbours[0][n], neighbours[1][n]);

        // simulate fire
        float prob = spread_onepix_prob_cpp(
          veg_target,
          terrain_burning,
          terrain_neighbour,
          coef_veg,
          coef_terrain,
          n,           // position argument, from 0 to 7
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

// a similar function only returning a burned_layer with an integer by burn 
// cycle, used to animate the spread.

// [[Rcpp::export]]
IntegerMatrix simulate_fire_animate(
    IntegerMatrix vegetation,
    arma::fcube terrain,
    IntegerMatrix ignition_cells,
    arma::frowvec coef,
    int n_veg_types = 6,
    float upper_limit = 1.0
) {
  
  // separate vegetation and terrain coefficients
  arma::frowvec coef_veg(n_veg_types);
  arma::frowvec coef_terrain(4);
  coef_veg = coef.subvec(0, n_veg_types - 1);
  coef_terrain = coef.subvec(n_veg_types, coef.size() - 1);
  
  // define landscape dimensions
  int n_row = vegetation.nrow();
  int n_col = vegetation.ncol();
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
  
  // The burned_step matrix will indicate the step at which each pixel was
  // burned, starting from the ignition (there will always be at least a "1"
  // pixel). Pixels with 0 are unburned.
  IntegerMatrix burned_step(n_row, n_col);
  
  // initialize with ignition_cells
  for(int i = 0; i <= end; i++) {
    burned_step(ignition_cells(0, i), ignition_cells(1, i)) = 1;
  }
  
  // initialize burning step
  int step = 2;
  
  while(burning_size > 0) {
    // Loop over all the burning cells to burn their neighbours. Use end_forward
    // to update the last position in burned_ids within this loop, without
    // compromising the loop's integrity.
    int end_forward = end;
    
    // Loop over burning cells in the cycle
    
    // b is going to keep the position in burned_ids that have to be evaluated
    // in this burn cycle
    for(int b = start; b <= end; b++) {
      
      // Get burning_cell's data
      arma::frowvec terrain_burning = terrain.tube(burned_ids(0, b), burned_ids(1, b));
      
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
        
        // Get vegetation class to know whether the cell is burnable
        int veg_target = vegetation(neighbours[0][n], neighbours[1][n]);
        
        // Is the cell burnable?
        bool burnable_cell = 
          (burned_step(neighbours[0][n], neighbours[1][n]) == 0) & // not burned
          (veg_target < 99);                                      // burnable
        if(!burnable_cell) continue;
        
        // obtain data from the neighbour
        arma::frowvec terrain_neighbour = terrain.tube(neighbours[0][n], neighbours[1][n]);
        
        // simulate fire
        float prob = spread_onepix_prob_cpp(
          veg_target,
          terrain_burning,
          terrain_neighbour,
          coef_veg,
          coef_terrain,
          n,           // position argument, from 0 to 7
          upper_limit
        );
        
        int burn = R::rbinom(1.0, prob);
        if(burn == 0) continue;
        
        // If burned,
        // store id of recently burned cell and
        // set step in burned_step
        // (but advance end_forward first)
        end_forward += 1;
        burned_ids(0, end_forward) = neighbours[0][n];
        burned_ids(1, end_forward) = neighbours[1][n];
        burned_step(neighbours[0][n], neighbours[1][n]) = step;
        
      } // end loop over neighbours of burning cell b
      
    } // end loop over burning cells from this cycle
    
    // update start and end
    start = end + 1;
    end = end_forward;
    burning_size = end - start + 1;
    
    // update step
    step += 1;
    
  } // end while
  
  return burned_step;
} 

// -----------------------------------------------------------------------

// The same function to be exported to R, only returning the binary burned_bin
// matrix.
// [[Rcpp::export]]
IntegerMatrix simulate_fire_cpp(
    IntegerMatrix vegetation,
    arma::fcube terrain,
    IntegerMatrix ignition_cells,
    arma::frowvec coef,
    int n_veg_types = 6,
    float upper_limit = 1.0
) {
  return simulate_fire_internal(
    vegetation,
    terrain,
    ignition_cells,
    coef,
    n_veg_types,
    upper_limit
  ).burned_bin;
}

// -----------------------------------------------------------------------

// The same function but deterministic, to test if the discrepancy between R and
// cpp is caused by seed problems

// [[Rcpp::export]]
IntegerMatrix simulate_fire_deterministic_cpp(
    IntegerMatrix vegetation,
    arma::fcube terrain,
    IntegerMatrix ignition_cells,
    arma::frowvec coef,
    int n_veg_types = 6,
    float upper_limit = 1.0
  ) {
  
  return simulate_fire_internal(
    vegetation,
    terrain,
    ignition_cells,
    coef,
    n_veg_types,
    upper_limit,
    [](double _, double x) { return (double)(x >= 0.5); }
  ).burned_bin;
}

// -------------------------------------------------------------------------

// simulate_fire_compare: same as simulate_fire, but returning objects to 
// compute discrepancy or similarity metrics (as a List):
//   IntegerMatrix burned_layer: binary matrix indicating the burned pixels,
//   IntegerMatrix burned_ids: ids as [row, col] column-vectors indicating the
//     position of burned pixels,
//   NumericVector counts_veg: burned pixels by veg_type. It has to be numeric
//     and not integer to allow non-integer divisions later.

typedef struct _s_burned_compare {
  IntegerMatrix burned_layer;
  IntegerMatrix burned_ids;
  NumericVector counts_veg; // need to be numeric to compute divisions later
} burned_compare;

burned_compare simulate_fire_compare_cpp(
    IntegerMatrix vegetation,
    arma::fcube terrain,
    IntegerMatrix ignition_cells,
    arma::frowvec coef,
    int n_veg_types = 6,
    float upper_limit = 1.0
) {

  burned_res burned = simulate_fire_internal(
    vegetation,
    terrain,
    ignition_cells,
    coef,
    n_veg_types,
    upper_limit
  );

  IntegerMatrix burned_bin = burned.burned_bin;
  IntegerMatrix burned_ids = burned.burned_ids;
  int end = burned.end;

  // Compute burned area by vegetation type
  NumericVector counts_veg(n_veg_types);
  for(int i = 0; i <= end; i++) {
    int veg_i = vegetation(burned_ids(0, i), burned_ids(1, i));
    counts_veg[veg_i] += 1;
  }
  
  return {burned_bin, burned_ids(_, seq(0, end)), counts_veg};
}

// [[Rcpp::export]]
List simulate_fire_compare(
    IntegerMatrix vegetation,
    arma::fcube terrain,
    IntegerMatrix ignition_cells,
    arma::frowvec coef,
    int n_veg_types = 6,
    float upper_limit = 1.0
) {
  
  burned_compare burned_com = simulate_fire_compare_cpp(
    vegetation,
    terrain,
    ignition_cells,
    coef,
    n_veg_types,
    upper_limit
  );

  // List to return:
  List L = List::create(Named("burned_layer") = burned_com.burned_layer,
                        Named("burned_ids") = burned_com.burned_ids,
                        Named("counts_veg") = burned_com.counts_veg);

  return L;
}