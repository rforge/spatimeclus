#include "STCXEM.h"

//[[Rcpp::export]]
S4  SpaTimeClusCpp(S4 input, List inputparam, NumericMatrix matT){
  cout << " deb c++" << endl;
  S4 * input_p=&input;
  STCXEM xem(input, inputparam, matT);
  xem.Run();
  //xem.Output(input_p);
  cout << " fin c++" << endl;
  return input;
}   