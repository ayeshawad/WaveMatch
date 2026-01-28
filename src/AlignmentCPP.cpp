#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix createAlignmentMatrix(int nrow, int ncol) {
  /*Filler*/
  NumericMatrix mat(nrow, ncol);
  for(int i = 0; i < nrow; i++) {
    for(int j = 0; j < ncol; j++) {
      mat(i, j) = 0;
    }
  }
  return mat;
}