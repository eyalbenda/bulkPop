#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector portMap(IntegerVector xpos, IntegerVector ypos, NumericVector xmap) {
  NumericVector ymap(ypos.size());
  ymap[0] = xmap[0];
  if(sum(ypos-ypos.sort())!=0)
  {
    stop("variants in Y are not sorted based on position!");
  }
  int curX = 0;
  for(int i=1;i<ypos.size();i++)
  {
    while(ypos[i]>xpos[curX] & xpos[curX]<max(xpos))
    {
      curX++;
    }
    ymap[i] = xmap[curX];
    Rcpp::checkUserInterrupt();
  }
  return ymap;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
//
// /*** R
// timesTwo(42)
// */
