#ifndef STCXEM_H
#define STCXEM_H


#include "STCdata.h"
#include "STCmodel.h"
#include "STCtune.h"
#include "STCparam.h"

class STCXEM{
  public:
  STCdata * m_data_p ;
  STCmodel * m_model_p ;
  STCtune * m_tune_p ;
  STCparam * m_paramCurrent_p ; 
  
  vector<STCparam> m_paramlist ;
  Col<double> m_loglikeSmall ;
  
  Mat<double> m_matT, m_tig, m_hessian;
  Col<double> m_poidspolynom;
  Cube<double>  m_Mjte;
  vector< vector < Cube<double> > >  m_sig ;
  
  bool m_nondegeneracy;
  
  STCXEM(){};
  STCXEM(const S4 &, const List &, const NumericMatrix &);
  ~STCXEM(){};  
  
  void Run();
  void SwitchParamCurrent(int);
  double ComputeLogLike();
  void Estep();
  void Mstep();
  void OneEM(const int, const double);
  void NewtonLogitWeighted(const int);
  void Output(S4 *);
};
#endif