#ifndef STCdata_H
#define STCdata_H

#include <iostream>
#include <iomanip>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class STCdata{
  public:
    Mat<double> m_x, m_map;
    int m_TT, m_JJ, m_n;
  
  STCdata(){};
  STCdata(const S4 & obj){
    this->m_x = as<mat>(obj.slot("x"));
    this->m_map  = as<mat>(obj.slot("map")); 
    this->m_TT = obj.slot("TT");
    this->m_JJ = obj.slot("JJ");
    this->m_n = obj.slot("n");
  };
  ~STCdata(){};  
};
#endif