#include <Rcpp.h>
using namespace Rcpp;

/*
 * Discrepancy functions to compare fires, mainly observed and simulated ones.
 *
 * Fires are represented as binary vectors [0, 1], taking 1 the burned pixels
 * and 0 the unburned ones.
 *
 * The main discrepancy functions are
 *  overlap
 *  delta_m
 *
 * Then, a wrapper function is used to compute the discrepancy considering many
 * fires (with a reference one or pairwise comparisons).
 *
 *  compare_fires
 */


// -------------------------------------------------------------------------


//' @title overlap_sp
//' @description Computes the overlap_sp index between two fires (binary vectors).
//'   (This is the ---- metric). fire_a and fire_b must have the same length and take
//'   only {0, 1}.
//' @return double: overlap_sp index, a value between 0 and 1, taking 1 when fires
//'   area identical and 0 when they are completely disjoint.
//' @param IntegerVector fire_a: a fire, binary vector with zeroes in unburned
//'   pixels and ones in burned pixels.
//' @param IntegerVector fire_b: another fire.

// [[Rcpp::export]]
double overlap_sp(IntegerVector fire_a, IntegerVector fire_b) {

  // check both fires have burned pixels
  if((sum(fire_a) == 0) | (sum(fire_b) == 0)) {
    stop("At least one fire has zero burned pixels.");
  }

  // number of cells burned in both fires
  int common = 0;
  for(int i = 0; i < fire_a.length(); i++) {
    common += ((fire_a[i] + fire_b[i]) == 2);
  }
  // burned cells in layers a and b
  int burned_a = sum(fire_a);
  int burned_b = sum(fire_b);

  // spatial overlap index
  double overlap_index = common / (burned_a + burned_b - common);

  return overlap_index;
}

/*** R
# test overlap_sp

# must be 1
fire_1 <- fire_2 <- rep(1, 6)
overlap_sp(fire_1, fire_2)

# must be 0:
fire_1 <- rep(c(1, 0), each = 3)
fire_2 <- rep(c(0, 1), each = 3)
overlap_sp(fire_1, fire_2)
*/


// -------------------------------------------------------------------------


//' @title delta_m
//' @description Computes the discrepancy measure proposed by Morales et al.
//'   (2015) between two fires (binary vectors). This measure is not symmetric,
//'   so a fire has to be set as the reference one.
//'   fire_a and fire_b must have the same length and take only {0, 1}.
//' @return double: discrepancy index.
//' @param IntegerVector fire_a: a fire, binary vector with zeroes in unburned
//'   pixels and ones in burned pixels. This is the reference fire (usually, the
//'   observed one in an ABC setting).
//' @param IntegerVector fire_b: another fire.
//' @param int veg_types: number of vegetation classes in the landscape.
//'   Defaults to 4.
//' @param NumericMatrix landscape: matrix with environmental data. This is used
//'   to get the burned area by vegetation type for each fire. The first
//'   veg_types - 1 columns must have binary vegetation layers, and one
//'   vegetation type is taken as background.
//' @param bool squared: whether to return the squared discrepancy or not.
//'   Defaults to false.

// (not used)
//' @param IntegerVector burnable: binary vector with ones in the burnable
//'   pixels in the landscape, used to get the shrubland layer.

// [[Rcpp::export]]
double delta_m(IntegerVector fire_a,
               IntegerVector fire_b,
               NumericMatrix landscape,
               int veg_types = 4,
               bool squared = 0) {

  // check both fires have burned pixels
  if((sum(fire_a) == 0) | (sum(fire_b) == 0)) {
    stop("At least one fire has zero burned pixels.");
  }

  // overall discrepancy (1 - overlap_spatial)
  double overall = 1 - overlap_sp(fire_a, fire_b);

  // get burned area by fire and vegetation type
  IntegerVector burned_veg_a(veg_types);
  IntegerVector burned_veg_b(veg_types);

  // Fill all vegetation types except the background one (shrubland)
  for(int v = 0; v < (veg_types - 1); v++) {
    int count_a = 0;
    int count_b = 0;
    for(int i = 0; i < fire_a.length(); i++) {
      if((landscape(i, v) == 1) & (fire_a[i] == 1)) count_a += 1;
      if((landscape(i, v) == 1) & (fire_b[i] == 1)) count_b += 1;
    }
    burned_veg_a[v] = count_a;
    burned_veg_b[v] = count_b;
  }
  // Fill background veg_type (shrubland)
  burned_veg_a[veg_types - 1] = sum(fire_a) - sum(burned_veg_a[seq(1, 3)]);
  burned_veg_b[veg_types - 1] = sum(fire_b) - sum(burned_veg_b[seq(1, 3)]);

  // compute delta.
  double delta = overall;
  for(int v = 0; v < veg_types; v++) {
    if(burned_veg_a[v] > 0) {
      delta += abs(burned_veg_a[v] - burned_veg_b[v]) / burned_veg_a[v];
    } else {
      delta += abs(burned_veg_a[v] - burned_veg_b[v]);
    }
  }

  return delta;
  if(squared == 1) return(pow(delta, 2.0));
}

/*** R
# test delta_m

# must be 0:
fire_1 <- fire_2 <- rep(1, 6)
landscape <- matrix(0, 6, 3) # all shrubland
delta_m(fire_1, fire_2, landscape)

# must be > 0:
fire_1 <- rep(c(1, 0), each = 3)
fire_2 <- rep(c(0, 1), each = 3)
delta_m(fire_1, fire_2, landscape)
*/


// -------------------------------------------------------------------------

//' @title overlap_vd
//' @description Computes the overlap between the vegetation type distribution
//'   ("vd") within each fire.
//' @param IntegerVector fire_a: a fire, binary vector with zeroes in unburned
//'   pixels and ones in burned pixels.
//' @param IntegerVector fire_b: another fire.
//' @param int veg_types: number of vegetation classes in the landscape.
//'   Defaults to 4.
//' @param NumericMatrix landscape: matrix with environmental data. This is used
//'   to get the burned area by vegetation type for each fire. The first
//'   veg_types - 1 columns must have binary vegetation layers, and one
//'   vegetation type is taken as background.

// [[Rcpp::export]]
double overlap_vd(IntegerVector fire_a,
                  IntegerVector fire_b,
                  NumericMatrix landscape,
                  int veg_types = 4) {

  // check both fires have burned pixels
  if((sum(fire_a) == 0) | (sum(fire_b) == 0)) {
    stop("At least one fire has zero burned pixels.");
  }

  // get burned area by fire and vegetation type
  IntegerVector burned_veg_a(veg_types);
  IntegerVector burned_veg_b(veg_types);

  // Fill all vegetation types except the background one (shrubland)
  for(int v = 0; v < (veg_types - 1); v++) {
    int count_a = 0;
    int count_b = 0;
    for(int i = 0; i < fire_a.length(); i++) {
      if((landscape(i, v) == 1) & (fire_a[i] == 1)) count_a += 1;
      if((landscape(i, v) == 1) & (fire_b[i] == 1)) count_b += 1;
    }
    burned_veg_a[v] = count_a;
    burned_veg_b[v] = count_b;
  }
  // Fill background veg_type (shrubland)
  burned_veg_a[veg_types - 1] = sum(fire_a) - sum(burned_veg_a[seq(1, 3)]);
  burned_veg_b[veg_types - 1] = sum(fire_b) - sum(burned_veg_b[seq(1, 3)]);

  // Get vegetation distribution by fire (normalized burned areas)
  NumericVector burned_dist_a(veg_types);
  NumericVector burned_dist_b(veg_types);

  for(int v = 0; v < veg_types; v++) {
    burned_dist_a[v] = burned_veg_a[v] / sum(burned_veg_a);
    burned_dist_b[v] = burned_veg_b[v] / sum(burned_veg_b);
  }

  // compute vegetation distribution overlap
  double overlap_index = 0;
  for(int v = 0; v < veg_types; v++) {
    overlap_index += std::min(burned_dist_a[v], burned_dist_b[v]);
  }

  return overlap_index;
}

/*** R
# test overlap_vd

# must be 1:
fire_1 <- fire_2 <- rep(1, 6)
landscape <- matrix(0, 6, 3) # all shrubland
overlap_vd(fire_1, fire_2, landscape)

# must be 1:
fire_1 <- rep(c(1, 0), each = 3)
fire_2 <- rep(c(0, 1), each = 3)
overlap_vd(fire_1, fire_2, landscape)

# must be 0:
landscape2 <- landscape
landscape2[1:3, 1] <- 1 # three pixels are subalpine
overlap_vd(fire_1, fire_2, landscape2) # must be 1
*/

// -------------------------------------------------------------------------

//' @title overlap_spvd
//' @description Computes a similarity measure defined as a weighted average
//'   between the spatial overlap ("sp") and the vegetation distribution overlap
//'   ("vd"). fire_a and fire_b must have the same length and take only {0, 1}.
//' @return double: similarity or discrepancy index.
//' @param IntegerVector fire_a: a fire, binary vector with zeroes in unburned
//'   pixels and ones in burned pixels.
//' @param IntegerVector fire_b: another fire.
//' @param int veg_types: number of vegetation classes in the landscape.
//'   Defaults to 4.
//' @param NumericMatrix landscape: matrix with environmental data. This is used
//'   to get the burned area by vegetation type for each fire. The first
//'   veg_types - 1 columns must have binary vegetation layers, and one
//'   vegetation type is taken as background.
//' @param NumericVector weights: non-negative weights for spatial overlap and
//'   vegetation distribution overlap, respectively. May not be normalized.
//' @param bool discrepancy: if true, a discrepancy, and not similarity measure
//'   is reported (1 - overlap).
//' @param bool squared: whether to return the squared discrepancy or not.
//'   Defaults to false. Ignored if discrepancy = FALSE.

// [[Rcpp::export]]
double overlap_spvd(IntegerVector fire_a,
                    IntegerVector fire_b,
                    NumericMatrix landscape,
                    //NumericVector weights = {0.8, 0.2},
                    int veg_types = 4,
                    bool discrepancy = 0,
                    bool squared = 0) {

  NumericVector weights = {0.8, 0.2};

  // get overlaps
  double overlap_spatial = overlap_sp(fire_a, fire_b);
  double overlap_vegdist = overlap_vd(fire_a, fire_b, landscape, veg_types);

  // get overall overlap (weighted average)
  double overlap_index = weights[0] / sum(weights) * overlap_spatial +
                         weights[1] / sum(weights) * overlap_vegdist;

  return overlap_index;

  if (discrepancy == 1) {
    if(squared == 1) {
      return(pow(1 - overlap_index, 2.0));
    } else {
      return(1 - overlap_index);
    }
  }

}
