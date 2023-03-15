#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include <spread_functions.cpp>

typedef struct _s_compare_result {
  double overlap_sp;

  double overlap_vd;
  double overlap_norm;
  double overlap_expquad;
  double overlap_quad;

  double sp_norm_5050;
  double sp_norm_7525;
  double sp_expquad_5050;
  double sp_expquad_7525;
  double sp_quad_5050;
  double sp_quad_7525;
} compare_result;

compare_result compare_fires_try_cpp(
    burned_compare fire1, burned_compare fire2,
    double lscale = 0.2) {

  // Extract list elements ------------------------------------------------

  IntegerMatrix burned1 = fire1.burned_layer;
  IntegerMatrix burned2 = fire2.burned_layer;

  IntegerMatrix burned_ids1 = fire1.burned_ids;
  IntegerMatrix burned_ids2 = fire2.burned_ids;

  double size1 = burned_ids1.ncol();
  double size2 = burned_ids2.ncol();

  NumericVector counts1 = fire1.counts_veg;
  NumericVector counts2 = fire2.counts_veg;

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

  compare_result indexes = {
    // pure indices
    .overlap_sp      = overlap_sp,

    .overlap_vd      = overlap_vd,
    .overlap_norm    = overlap_norm,
    .overlap_expquad = overlap_expquad,
    .overlap_quad    = overlap_quad,

    // mixture indices
    .sp_norm_5050    = (0.50 * overlap_sp + 0.50 * overlap_norm),
    .sp_norm_7525    = (0.75 * overlap_sp + 0.25 * overlap_norm),
    .sp_expquad_5050 = (0.50 * overlap_sp + 0.50 * overlap_expquad),
    .sp_expquad_7525 = (0.75 * overlap_sp + 0.25 * overlap_expquad),
    .sp_quad_5050    = (0.50 * overlap_sp + 0.50 * overlap_quad),
    .sp_quad_7525    = (0.75 * overlap_sp + 0.25 * overlap_quad)
  };

  return indexes;
}

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

  // Extract list elements ------------------------------------------------

  IntegerMatrix burned1 = fire1["burned_layer"];
  IntegerMatrix burned2 = fire2["burned_layer"];

  IntegerMatrix burned_ids1 = fire1["burned_ids"];
  IntegerMatrix burned_ids2 = fire2["burned_ids"];

  NumericVector counts1 = fire1["counts_veg"];
  NumericVector counts2 = fire2["counts_veg"];

  compare_result indexes = compare_fires_try_cpp(
    {burned1, burned_ids1, counts1},
    {burned2, burned_ids2, counts2}
  );

  return NumericVector::create(
    // pure indices
    Named("overlap_sp")      = indexes.overlap_sp,

    Named("overlap_vd")      = indexes.overlap_vd,
    Named("overlap_norm")    = indexes.overlap_norm,
    Named("overlap_expquad") = indexes.overlap_expquad,
    Named("overlap_quad")    = indexes.overlap_quad,

    // mixture indices
    Named("sp_norm_5050")    = indexes.sp_norm_5050,
    Named("sp_norm_7525")    = indexes.sp_norm_7525,
    Named("sp_expquad_5050") = indexes.sp_expquad_5050,
    Named("sp_expquad_7525") = indexes.sp_expquad_7525,
    Named("sp_quad_5050")    = indexes.sp_quad_5050,
    Named("sp_quad_7525")    = indexes.sp_quad_7525
  );
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
//'   metric (columns) by simulated fire (rows)

//' @param arma::cube landscape: Environmental data from the whole landscape.
//'   See description in spread_around. The 3rd dimension contains each layer
//'   (covariates). The landscape is represented in matrix form.
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

// [[Rcpp::export]]
NumericMatrix emulate_loglik_try(
    arma::cube landscape,
    IntegerMatrix ignition_cells,
    IntegerMatrix burnable,
    arma::rowvec coef,
    int wind_layer,
    int elev_layer,
    arma::rowvec distances,
    double upper_limit,
    List fire_ref,
    int n_replicates = 10
) {

  CharacterVector names = CharacterVector::create(
    "overlap_sp",
    "overlap_vd",
    "overlap_norm",
    "overlap_expquad",
    "overlap_quad",

    "sp_norm_5050",
    "sp_norm_7525",
    "sp_expquad_5050",
    "sp_expquad_7525",
    "sp_quad_5050",
    "sp_quad_7525"
  );

  NumericMatrix similarity(n_replicates, names.length());

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

  colnames(similarity) = names;
  return similarity;
}


// The same function but evaluating many particles all at once. It's implemented
// this way because turning R objects into C++ ones takes a long time, so, if a 
// landscape is going to be used to simulate many fires from many particles, 
// translating the data only once should be more efficient. 

// The coef argument is now replaced by an arma::mat, with particles in rows and
// parameters in columns.

//' @param arma::mat particles: parameters row_vectors (particles) in logistic 
//' regression to compute the spread probability as a function of covariates.

// [[Rcpp::export]]
arma::cube emulate_loglik_try_par(
    arma::cube landscape,
    IntegerMatrix ignition_cells,
    IntegerMatrix burnable,
    arma::mat particles, // it was coef before, only a vector, now a matrix
    int wind_layer,
    int elev_layer,
    arma::rowvec distances,
    double upper_limit,
    List fire_ref,
    int n_replicates = 10
) {
  
  int n_particles = particles.n_rows;
  
  CharacterVector names = CharacterVector::create(
    "overlap_sp",
    "overlap_vd",
    "overlap_norm",
    "overlap_expquad",
    "overlap_quad",
    
    "sp_norm_5050",
    "sp_norm_7525",
    "sp_expquad_5050",
    "sp_expquad_7525",
    "sp_quad_5050",
    "sp_quad_7525"
  );
  
  arma::cube similarity(n_replicates, names.length(), n_particles);
  arma::rowvec ind_temp(names.length()); 
  
  for(int part = 0; part < n_particles; part++) {
    for(int i = 0; i < n_replicates; i++) {
      
      // simulate_fire
      List fire_sim = simulate_fire_compare(
        landscape,
        ignition_cells,
        burnable,
        particles.row(part), // this is coef argument, if coef = ..., it doesint owrk
        wind_layer,
        elev_layer,
        distances,
        upper_limit
      );
      
      ind_temp = compare_fires_try(fire_ref, fire_sim);
      for(int j = 0; j < names.length(); j++) {
        similarity(i, j, part) = ind_temp[j];
      } 
      // filling the array in vectorized way did not work: 
      
      // similarity.slice(part).row(i).fill(ind_temp);
      // similarity.slice(part).row(i).fill(compare_fires_try(fire_ref, fire_sim));
    }
  }
  
  // colnames(similarity) = names; // no anda
  return similarity;
}