#include "STCXEM.h"

//[[Rcpp::export]]
S4  SpaTimeClusCpp(S4 input, List inputparam, NumericMatrix matT, NumericMatrix toollogistic){
  S4 * input_p=&input;
  STCXEM xem(input, inputparam, matT, toollogistic);
  xem.Run();
  xem.Output(input_p);
  return input;
}   