#include "STCXEM.h"

//[[Rcpp::export]]
S4  SpaTimeClusCpp(S4 input, List inputparam, NumericMatrix matT){
  S4 * input_p=&input;
  STCXEM xem(input, inputparam, matT);
  xem.Run();
  xem.Output(input_p);
  return input;
}   