#include <Rcpp.h>
using namespace Rcpp;

/*
 * Sample particles without replacement using weights, used for the SMC move with
 * DE-MCMC. It requires 3 distinct particles, but sampled with weights.
 */

// [[Rcpp::export]]
IntegerMatrix sample_triplets_weighted(const NumericVector& w) {
  int N = w.size();
  IntegerMatrix out(N, 3);
  
  NumericVector prob(N);

  for (int i = 0; i < N; ++i) {
    
    // copy weights
    std::copy(w.begin(), w.end(), prob.begin());
    
    double total_w = std::accumulate(prob.begin(), prob.end(), 0.0);
    
    for (int k = 0; k < 3; ++k) {
      
      double u = R::runif(0.0, total_w);
      double cumsum = 0.0;
      int chosen = -1;
      
      for (int j = 0; j < N; ++j) {
        if (prob[j] <= 0.0) continue;
        cumsum += prob[j];
        if (u <= cumsum) {
          chosen = j;
          break;
        }
      }
      
      // store 1-based index
      out(i, k) = chosen + 1;
      
      // remove chosen element
      total_w -= prob[chosen];
      prob[chosen] = 0.0;
    }
  }
  
  return out;
}
