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
 *
 *
 *  -------------------------------------------------------------
 *
 *  In a few functions I define double or NumericVector variables that could be
 *  integer, but I do it so for computations involving these integers to return
 *  double results. Apparently, double x = int1 + int1 is not double, but int.
 *
 *  -------------------------------------------------------------
 *  
 *  
 *  
 *                      IMPORTANT  !!!!!!!
 *  
 *  MASSIVE IMPROVEMENTS TO DO, BASED ON THE MATRIX REPRESENTATION OF THE 
 *  LANDSCAPE. ALL FUNCTIONS BELOW EVALUATE THE WHOLE LANDSCAPE 
 *  UNNECESSARILY.
 *  
 *  Improvements:
 *  
 *  With the matrix representation of the landscape, where burning (burned)
 *  cell ids are stored in a matrix, using the burning_ids we can avoid evaluating
 *  the whole landscape. 
 *  
 *  To count burned by veg_type we could do this:
 *  
 *  // end is the last position filled in burning_ids, so its the number of 
 *  // pixels burned. Check if it should be summed one or not to be the size.
 *  
 *  NumericVector counts(veg_types);
 *  
 *  for(int i = 0; i <= end; i++) {
 *    for(int v = 1; v < veg_types; v++) { // starts from 1 because 0 is for shrubland
 *      counts(v) += landscape(burning_ids(0, i), burning_ids(1, i), v-1);
 *    }
 *    // landscape must have the vegetation layers in the first veg_types - 1 layers 
 *  }
 *  
 *  // fill shrubland
 *  counts(v) = (end + 1) - sum(counts(seq(1, veg_types - 1)));
 *  
 *  NumericVector props(veg_types);
 *  for(int i = 0; i < veg_types; i++) props(i) = counts(i) / end;
 *  
 *  
 *  
 *  // ---------------------------
 *  
 *  Another improvement for overlap_sp would be to evaluate the intersection
 *  using only the smaller fire. So, only in those ids the intersection counter
 *  would operate:
 *  
 *  int smaller_size = min({end, reference_fire_size});
 *  
 *  double both = 0.0;
 *  
 *  for(int i = 0; i < smaller_size; i++) {
 *    both(i) += burned_layer_larger(burning_ids_smaller(0, i), 
 *                                   burning_ids_smaller(1, i));
 *  }
 *  
 *  // this subsetting could probably be improved using armadillo arrays
 *  
 *  overlap = both / (size_a + size_b - both)
 *  
 *  
 *  //// COMPARE WHETHER A COUNTER LIKE THOSE ABOVE IS SLOWER THAN 
 *  sum(vector[subset]), to test if using armadillo structures is better.
 *  
 *  
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
  double common = 0;                                 // changed int by double
  for(int i = 0; i < fire_a.length(); i++) {
    common += ((fire_a[i] + fire_b[i]) == 2);
  }
  // burned cells in layers a and b
  double burned_a = sum(fire_a);                                 // changed int by double
  double burned_b = sum(fire_b);                                 // changed int by double
  
  // spatial overlap index
  double overlap_index = common / (burned_a + burned_b - common);
  /*
   * IMPORTANT: in this computation all variables must be double for the result
   * to be double. If not, the result is integer, although it's defined as double.
   */
  
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

# must be between 0 and 1:
fire_1 <- rep(c(1, 0), each = 3)
fire_2 <- c(rep(1, 5), 0)
overlap_sp(fire_1, fire_2)

*/


// -------------------------------------------------------------------------


//' @title delta_m
//' @description Computes the discrepancy measure proposed by Morales et al.
//'   (2015) between two fires (binary vectors). Morales used the area of
//'   fire_a in the denominator, but to make the distance symmetric, it
//'   is divided by mean(area(fire_a), area(fire_b)).
//'   fire_a and fire_b must have the same length and take only {0, 1}.
//'   delta_m_raw is almost the same, but does not add the complement of
//'   overlap. This allows to know how much does each component weight in
//'   delta_m.
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
    
    double mean_area = (burned_veg_a[v] + burned_veg_b[v]) / 2.0;
    
    if(mean_area > 0) {
      delta += abs(burned_veg_a[v] - burned_veg_b[v]) / mean_area;
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


// [[Rcpp::export]]
double delta_m_raw(IntegerVector fire_a,
                   IntegerVector fire_b,
                   NumericMatrix landscape,
                   int veg_types = 4,
                   bool squared = 0) {
  
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
  
  // compute delta.
  double delta = 0;
  for(int v = 0; v < veg_types; v++) {
    
    double mean_area = (burned_veg_a[v] + burned_veg_b[v]) / 2.0;
    
    if(mean_area > 0) {
      delta += abs(burned_veg_a[v] - burned_veg_b[v]) / mean_area;
    }
  }
  
  return delta;
  if(squared == 1) return(pow(delta, 2.0));
}



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
  NumericVector burned_veg_a(veg_types);  // changed integer to numeric class
  NumericVector burned_veg_b(veg_types);
  
  // Fill all vegetation types except the background one (shrubland)
  for(int v = 0; v < (veg_types - 1); v++) {
    double count_a = 0;
    double count_b = 0;
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
  
  /*
   * IMPORTANT: in this computation all variables must be double for the result
   * to be double. If not, the result is integer, although it's defined as double.
   */
  
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
overlap_vd(fire_1, fire_2, landscape2)

# must be between 0 and 1:
fire_1 <- rep(c(1, 0), each = 3)
fire_2 <- c(1, 1, 1, 1, 0, 0)
overlap_vd(fire_1, fire_2, landscape2)

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

/*** R
# # test overlap_spvd
#
# # must be 1
# fire_1 <- fire_2 <- rep(1, 6)
# landscape <- matrix(0, 6, 3) # all shrubland
# overlap_spvd(fire_1, fire_2, landscape)
#
# # must be 0.3:
# fire_1 <- rep(c(1, 0), each = 3)
# fire_2 <- rep(c(0, 1), each = 3)
# overlap_spvd(fire_1, fire_2, landscape, weights = c(0.7, 0.3))
#
# # must be 0:
# landscape2 <- landscape
# landscape2[1:3, 1] <- 1 # three pixels are subalpine
# overlap_spvd(fire_1, fire_2, landscape2)
*/


// -------------------------------------------------------------------------

/*
 
 función para comparar muchos fuegos. hacerlo en R, estoy tardando too much
 
 //' @title compare_fires
 //' @description Computes similarity or dissimilarity measures for a matrix of
 //'   fires (each column is a fire, each row is a pixel) from the same landscape.
 //'   It allows to compare all fires pairwise or take one as the reference. In
 //'   that case, the first column is considered, and reference should be set to
 //'   TRUE.
 //' @return NumericVector: similarity or discrepancy indexes.
 //'   If reference = FALSE, the result is a diagonal matrix.
 //' @param IntegerMatrix fires: matrix with fires in columns, pixels in rows,
 //'   occurring on the same landscape.
 //' @param char metric: the similarity/dissimilarity measure desired, options are
 //'   ("overlap_sp", "overlap_vd", "overlap_spvd", "delta_m"). All except
 //'   the first require a landscape and veg_types specification.
 //' @param int veg_types: number of vegetation classes in the landscape.
 //'   Defaults to 4. (Only used for a few metrics.)
 //' @param NumericMatrix landscape: matrix with environmental data. This is used
 //'   to get the burned area by vegetation type for each fire. The first
 //'   veg_types - 1 columns must have binary vegetation layers, and one
 //'   vegetation type is taken as background. (Only used for a few metrics.)
 //' @param NumericVector weights: non-negative weights for spatial overlap and
 //'   vegetation distribution overlap, respectively. May not be normalized.
 //' @param bool discrepancy: if true, a discrepancy, and not similarity measure
 //'   is reported (1 - overlap). (Only used for a few metrics.)
 //' @param bool squared: whether to return the squared discrepancy or not.
 //'   Defaults to false. Ignored if discrepancy = FALSE.
 
 double compare_fires(IntegerMatrix fires,
 NumericMatrix landscape,
 NumericVector weights = {0.8, 0.2},
 String metric("overlap_sp"),
 bool reference = 0,
 int veg_types = 4,
 bool discrepancy = 0,
 bool squared = 0) {
 
 // compare fires when there are only two
 if(fires.ncol() < 3) {
 
 
 
 }
 
 }
 
 ///// ALGO NO ANDUVO Y LO ABANDONÉ
 
 */


// ---------------------------------------------------------------------

// Many discrepancy functions

/*
 * This takes as argument two lists that bear the burned_layer, 
 * burned_ids (0 indexing) and
 * counts_veg, the burned pixels by veg_type (shrubland, subalpine, wet, dry).
 */

// [[Rcpp::export]]
NumericVector compare_fires(List fire1, List fire2) {
  
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
  
  // overlap_i -----------------------------------------------------------
  
  // overlap index consdering difference in burned area by vegetation type, 
  // with i from iván
  double diff_prop_raw = 0.0;
  
  for(int v = 0; v < veg_types; v++) {
    double max_area = std::max(counts1[v], counts2[v]);
    if(max_area > 0.0) {
      diff_prop_raw += abs(counts1(v) - counts2(v)) / max_area;
    }
  }
  
  double veg_num = counts1.length();
  double overlap_i = 1.0 - diff_prop_raw / veg_num;
  
  // delta_m_raw -----------------------------------------------------------
  
  // discrepancy index consdering difference in burned area by vegetation type, 
  // with m from Morales. This overlap is based on the discrepancy function
  // used in Morales et al. 2015.
  
  double delta_m_raw = 0.0;
  for(int v = 0; v < veg_types; v++) {
    
    if(counts1[v] > 0) {
      delta_m_raw += abs(counts1[v] - counts2[v]) / counts1[v]; 
    } // first fire is used as reference
    else {
      delta_m_raw += abs(counts1[v] - counts2[v]);
    }
  }
  
  // the same, but using the mean area in the denumerator to make it symmetric
  double delta_m_mean = 0.0;
  for(int v = 0; v < veg_types; v++) {
    
    double mean_area = (counts1[v] + counts2[v]) / 2.0;
    
    if(mean_area > 0) {
      delta_m_mean += abs(counts1[v] - counts2[v]) / mean_area;
    }
  }
  
  // ---------------------------------------------------------------------
  
  NumericVector indexes = NumericVector::create(
    Named("overlap_sp") = overlap_sp,
    Named("overlap_vd") = overlap_vd,
    Named("overlap_i") = overlap_i,
    Named("delta_m_raw") = delta_m_raw,
    Named("delta_m_mean") = delta_m_mean
  );
  
  return indexes;
}

