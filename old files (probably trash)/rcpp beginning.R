library(Rcpp)

# cpptxt <- "
# int fibonacci(const int x) {
#   if (x < 2) return(x);
#   return (fibonacci(x - 1)) + fibonacci(x - 2);
# }"
# 
# fibCpp <- cppFunction(cpptxt)
# 
# fibCpp(6)

vec_code <- "
NumericMatrix makeMat(int rows, int cols, int fill) {
  NumericMatrix m(rows, cols);
  
  int total = rows * cols;
  // fill matrix
  for(int i = 0; i < rows; i++) {
    for(int j = 0; j < cols; j++) {
      m(i, j) = fill;
    }
  }
  
  return m;
}  
"

makeMat <- cppFunction(vec_code)


makeMat(3, 5, 111)
