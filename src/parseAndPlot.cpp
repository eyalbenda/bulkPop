#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector summarizePop(List hapPop,int maxLength) {
  if(maxLength<1)
  {
    throw std::invalid_argument("length of haplotype must be larger than 0");
  }
  NumericVector counts(maxLength);
  NumericVector size(maxLength);
  List::iterator it;
  for(it = hapPop.begin();it!=hapPop.end();++it)
  {
    LogicalVector Hap = as<LogicalVector>(*it);
    for(int i=0;i<Hap.length();i++)
    {
      if(Hap[i])
      {
        counts[i] += 1;
      }
      size[i]+=1;
    }
  }
  return counts/size;
}

