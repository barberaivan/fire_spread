library(Rcpp)

cppFunction(depends = "RcppArmadillo", ' 
arma::vec give_part_of_cube(arma::cube x) {
  arma::vec z(3);
  z = x.tube(2, 2);
  return z;
}
')

(ar1 <- array(1, dim = c(3, 3, 3)))
(v <- give_part_of_cube(ar1))
