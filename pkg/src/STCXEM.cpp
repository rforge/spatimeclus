#include "STCXEM.h"

STCXEM::STCXEM(const S4 & input, const List & inputparam, const NumericMatrix & matT, const NumericMatrix & toollogistic){
  m_data_p = new STCdata( as<S4>(input.slot("data")) );
  m_model_p = new STCmodel( as<S4>(input.slot("model")) );
  m_tune_p = new STCtune( as<S4>(input.slot("tune")) );
  m_paramCurrent_p = new STCparam( as<S4>(inputparam[0]) );
  m_paramlist.resize(m_tune_p->m_nbinitSmall);
  for (int ini=0; ini< m_tune_p->m_nbinitSmall ; ini++) m_paramlist[ini] = STCparam( as<S4>(inputparam[ini]) );
  m_loglikeSmall= ones<vec>(m_tune_p->m_nbinitSmall) * log(0);
  m_matT = as<mat>(matT);
  m_toollogistic = as<mat>(toollogistic);
  m_weightlogistic=cube(m_data_p->m_JJ * m_data_p->m_TT, m_model_p->m_K, m_model_p->m_G, fill::zeros);
  m_condintramargin=cube(m_data_p->m_JJ * m_data_p->m_TT, m_data_p->m_n, m_model_p->m_G, fill::zeros);
  m_sig.resize(m_model_p->m_G);
  for (int g=0; g<m_model_p->m_G; g++)  m_sig[g] = cube(m_data_p->m_JJ * m_data_p->m_TT, m_data_p->m_n, m_model_p->m_K, fill::zeros);
  m_tig= mat(m_data_p->m_n, m_model_p->m_G, fill::zeros);
  m_poidspolynom= mat(m_data_p->m_JJ * m_data_p->m_TT, m_model_p->m_K, fill::zeros);
  m_hessian=mat(4*(m_model_p->m_K-1), 4*(m_model_p->m_K-1), fill::zeros);
  m_degeneracy = 0;
}

void STCXEM::SwitchParamCurrent(int ini){m_paramCurrent_p = &m_paramlist[ini];}

void STCXEM::Run(){  
  for (int ini=0; ini< m_tune_p->m_nbinitSmall ; ini++){
    m_degeneracy=0;
    SwitchParamCurrent(ini);
    OneEM(m_tune_p->m_nbiterSmall, m_tune_p->m_tol);
    if (m_degeneracy){
      m_loglikeSmall(ini) = -99999999999999999;
    }else{
      m_loglikeSmall(ini) = ComputeLogLike();
    }    
  }
  uvec indices = sort_index(m_loglikeSmall);
  for (int tmp1=0; tmp1< m_tune_p->m_nbinitKept; tmp1++){
    m_degeneracy=0;
    SwitchParamCurrent(indices(m_tune_p->m_nbinitSmall - tmp1 - 1));
    OneEM(m_tune_p->m_nbiterKept, m_tune_p->m_tol);
    if (m_degeneracy){
      m_loglikeSmall(indices(m_tune_p->m_nbinitSmall - tmp1 - 1))  = -99999999999999999;
    }else{
     m_loglikeSmall(indices(m_tune_p->m_nbinitSmall - tmp1 - 1)) = ComputeLogLike();
    }  
  }
  m_degeneracy=0;
  uword  index;
  double indicebest = (m_loglikeSmall).max(index);
  SwitchParamCurrent(index);
  indices = sort_index(m_loglikeSmall);
}


double STCXEM::ComputeLogLike(){
  double mu=0, loglike=   -99999999999999999;
  if (m_degeneracy == 0){
    int repere=0;
    Mat<double> center=mat(m_data_p->m_TT, m_model_p->m_K, fill::zeros);
    Mat<double> tmpind=mat(m_data_p->m_JJ, m_data_p->m_n, fill::zeros);
    m_condintramargin=cube(m_data_p->m_JJ * m_data_p->m_TT, m_data_p->m_n, m_model_p->m_G, fill::zeros);
    m_tig.each_row() = trans(log(m_paramCurrent_p->m_proportions));
    for (int g=0; g<m_model_p->m_G;g++){
      m_poidspolynom = m_toollogistic * trans(m_paramCurrent_p->m_lambda[g]);
      m_poidspolynom = exp(m_poidspolynom.each_col() - max(m_poidspolynom,1));
      m_poidspolynom.each_col() /= sum(m_poidspolynom,1);
      center = m_matT * trans(m_paramCurrent_p->m_beta[g]);
      for (int k=0; k< m_model_p->m_K; k++){
        for (int tt=0; tt<m_data_p->m_TT; tt++){
          tmpind=exp(-0.5 * pow(m_data_p->m_x.submat(tt*m_data_p->m_JJ, 0, (tt+1)*m_data_p->m_JJ-1, m_data_p->m_n -1) -  center(tt, k), 2)/m_paramCurrent_p->m_sigma(g,k)) / sqrt( 2*M_PI*m_paramCurrent_p->m_sigma(g,k));
          m_sig[g].subcube(tt*m_data_p->m_JJ, 0, k, (tt+1)*m_data_p->m_JJ-1, m_data_p->m_n -1, k) =  tmpind.each_col() % m_poidspolynom.submat(tt*m_data_p->m_JJ, k, (tt+1)*m_data_p->m_JJ-1, k);
        } 
        m_condintramargin.slice(g) += m_sig[g].slice(k);
      }
      m_tig.col(g) += trans(sum(log(m_condintramargin.slice(g)),0)) ;
    }
    Col<double> maxiclasses = max(m_tig,1);
    loglike = sum(maxiclasses);
    m_tig = exp(m_tig.each_col() - maxiclasses);
    maxiclasses = sum(m_tig,1);
    loglike=loglike + sum(log(maxiclasses));
  }
  return loglike;
}


void STCXEM::Estep(){
  m_tig = m_tig.each_col() / sum(m_tig,1);
  for (int g=0; g<m_model_p->m_G; g++){
    for (int k=0; k< m_model_p->m_K; k++){
      m_sig[g].slice(k) = m_sig[g].slice(k) /  m_condintramargin.slice(g);
      m_sig[g].slice(k) = m_sig[g].slice(k).each_row() % trans(m_tig.col(g));
      if (any(m_tig.col(g) <= 0)) (m_sig[g].slice(k)).elem(find(m_condintramargin.slice(g) <= 0)).zeros();
      (m_weightlogistic.slice(g)).col(k) = sum(m_sig[g].slice(k),1);
    }
  }
}

void STCXEM::NewtonLogitWeighted(const int g){
  Col<double> u = sum(m_weightlogistic.slice(g), 1);  
  Mat<double> ratioexp = m_toollogistic * trans(m_paramCurrent_p->m_lambda[g]);
  ratioexp = exp(ratioexp.each_col() - max(ratioexp, 1));
  ratioexp.each_col() /= sum(ratioexp,1);
  for (int k1=0; k1< (m_model_p->m_K-1); k1++){
    for (int k2=k1; k2< (m_model_p->m_K-1); k2++){
      for (int h1=0; h1 < 4; h1++){
        for (int h2=0; h2 < 4; h2++){
          m_hessian(k1*4 + h1, k2*4 + h2) = sum(u % m_toollogistic.col(h1) % m_toollogistic.col(h2) % ratioexp.col(k1+1) % ratioexp.col(k2+1) );
          m_hessian(k2*4 + h2, k1*4 + h1) = m_hessian(k1*4 + h1, k2*4 + h2);
        }
      }
    }
    for (int h1=0; h1 < 4; h1++){
      for (int h2=h1; h2 < 4; h2++){
        m_hessian(k1*4 + h1, k1*4 + h2) -= sum(u % m_toollogistic.col(h1) % m_toollogistic.col(h2) % ratioexp.col(k1+1));
        m_hessian(k1*4 + h2, k1*4 + h1) = m_hessian(k1*4 + h1, k1*4 + h2);
      }
    }
  }
  // Maintenant u contient les valeurs de parametres updates en vecteur 
  Col<double> currentparam=trans(vectorise(m_paramCurrent_p->m_lambda[g].rows(1, m_model_p->m_K-1),1));
  Col<double> gradientcommun = vectorise(trans(m_toollogistic) * m_weightlogistic.slice(g).cols(1, m_model_p->m_K-1)) ;
  Col<double> gradientautre = vectorise(trans(u % m_toollogistic.each_col()) *  ratioexp.cols(1, m_model_p->m_K-1));
  bool okinversion = solve(u, m_hessian,  gradientcommun - gradientautre, solve_opts::no_approx);
  if (okinversion == 1){
    u = currentparam - u;
    for (int k=0; k< (m_model_p->m_K-1); k++) m_paramCurrent_p->m_lambda[g].row(k+1) = trans(u.subvec(k*4, (k*4)+3));
  }else{
    m_degeneracy = 1;
  }
}

void STCXEM::Mstep(){
  m_paramCurrent_p->m_proportions = trans(sum(m_tig,0) / m_data_p->m_n);
  for (int g=0; g<m_model_p->m_G; g++){
    for (int k=0; k< m_model_p->m_K; k++){
      Col<double> tmp=sum(m_sig[g].slice(k), 1);
      Col<double> tmpX = sum(m_sig[g].slice(k) % m_data_p->m_x, 1);
      Col<double> lambdamat = vec( m_data_p->m_TT, fill::zeros);
      Col<double> lambdaX= vec( m_data_p->m_TT, fill::zeros);
      for (int t=0; t< m_data_p->m_TT; t++){
        for (int j=0; j< m_data_p->m_JJ; j++){
          lambdamat(t) += tmp(m_data_p->m_JJ*t + j);
          lambdaX(t) += tmpX(m_data_p->m_JJ*t + j);
        }
      }
      Col <double> inv=trans(m_paramCurrent_p->m_beta[g].row(k));
      bool okinversion=solve(inv, trans(m_matT) * diagmat(lambdamat) * m_matT, trans(m_matT) * lambdaX, solve_opts::no_approx);
      if (okinversion == 1){
        m_paramCurrent_p->m_beta[g].row(k) = trans(inv);     
      }else{
        m_degeneracy = 1;
        goto stop;
      }
      Col<double> tmpmean= m_matT * trans(m_paramCurrent_p->m_beta[g].row(k));
      Col<double> mean= vec( m_data_p->m_TT*m_data_p->m_JJ, fill::zeros);
      for (int t=0; t< m_data_p->m_TT; t++) mean.subvec(m_data_p->m_JJ*t, m_data_p->m_JJ*(t+1)-1) = vec(m_data_p->m_JJ, fill::ones) * tmpmean(t);
      m_paramCurrent_p->m_sigma(g,k) = sum( sum(pow(m_data_p->m_x.each_col() - mean, 2) % m_sig[g].slice(k))) / sum( sum(m_sig[g].slice(k)) );
    }
    if (m_model_p->m_K>1) { NewtonLogitWeighted(g);}
  }
   stop:
   ;
  
}

void STCXEM::OneEM(const int itermax, const double tol){
  m_degeneracy=0;
  double loglike = ComputeLogLike(), prec = log(0);
  int it=0;
  while ( (it<itermax) && ((loglike-prec)>tol) ){
    it ++;
    Estep();
    Mstep();
    if (m_degeneracy==1){
      prec = log(0);
      loglike = prec;
    }else{
      prec = loglike;
      loglike = ComputeLogLike();
      if (loglike< (prec-tol)) cout << "Error in EM algorithm (loglikelihood decreases)" << endl << "diff: "<< loglike - prec<< " loglike: " << loglike<< " prec: " << prec << " tol" << tol << " iteration " << it << endl;
    }
  } 
}

void STCXEM::Output(S4 * reference_p){
  as<S4>(reference_p->slot("criteria")).slot("loglike") = wrap(ComputeLogLike());
  as<S4>(reference_p->slot("param")).slot("proportions") = wrap(trans(m_paramCurrent_p->m_proportions));
  as<S4>(reference_p->slot("param")).slot("lambda") = wrap(m_paramCurrent_p->m_lambda);
  as<S4>(reference_p->slot("param")).slot("sigma") = wrap(m_paramCurrent_p->m_sigma);
  as<S4>(reference_p->slot("param")).slot("beta") = wrap(m_paramCurrent_p->m_beta);
}