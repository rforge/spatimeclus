##' SpaTimeClus a package for clustering spatio-temporal data
##'
##' todo
##'
##' \tabular{ll}{
##'   Package: \tab SpaTimeClus\cr 
##'   Type: \tab Package\cr 
##'   Version: \tab 1.0.0\cr
##'   Date: \tab 2015-11-03\cr 
##'   License: \tab GPL-2\cr 
##'   LazyLoad: \tab yes\cr
##' }
##'
##' todo 
##' 
##' @name SpaTimeClus-package
##' @aliases SpaTimeClus
##' @rdname SpaTimeClus-package
##' @docType package
##' @keywords package
##' @import parallel
##' @import Rcpp
##' @import RcppArmadillo
##' @import methods
##' @importFrom stats runif
##' @useDynLib SpaTimeClus
##'
##' @author
##' Author: Marbac M., and McNicholas P.
##'
##' @references Marbac M., and McNicholas P. 
##' 
##' @examples
##' data(ech)
##' 
##' output <- spatimeclus(ech = ech, G = 2, K = 2, Q = 2, crit = "BIC", nbcores = 6)
##' 
##' summary(output)
##' 
##' print(output)
##' 
##' plot(output) 
##' 
NULL



spatimeclusModelKnown <- function(ech, model, tol=0.1, param=NULL, nbcores=1, nbinitSmall=100, nbinitKept=10, nbiterSmall=10, nbiterKept=500){
  nbcores <- min(detectCores(all.tests = FALSE, logical = FALSE),  nbcores)
  nbinitSmall <- nbcores * ceiling(nbinitSmall / nbcores) 
  nbinitKept <- nbcores * ceiling(nbinitKept / nbcores) 
  if (is.null(param)){
    param <- list()
    for (it in 1:nbinitSmall)  param[[it]] <- initparam(ech, model)
  }
  input <-new("STCresults", 
              model=model, 
              data=ech, 
              criteria=new("STCcriteria", loglike= -Inf),
              tune=new("STCtune", 
                       tol=tol, 
                       nbinitSmall=ceiling(nbinitSmall / nbcores), 
                       nbinitKept=ceiling(nbinitKept / nbcores), 
                       nbiterSmall=nbiterSmall, 
                       nbiterKept=nbiterKept
                       )
              )
  
  matT <- t(exp(sweep(log(matrix(1:ech@TT, model@Q+1, ech@TT, byrow=TRUE)), 1, 0:model@Q, "*")))
  toollogistic <- cbind(ech@map, rep(1, ech@JJ))
  for (tt in 2:ech@TT)    toollogistic <- rbind(toollogistic, cbind(ech@map, rep(tt,  ech@JJ)))
  toollogistic <- cbind(toollogistic, rep(1, nrow(toollogistic)))
  
  if (nbcores>1){
    paramparallel <- list()
    compteur <- 0
    for (c in 1:nbcores){
      paramparallel[[c]] <- list()
      for (i in 1:input@tune@nbinitSmall){
        compteur <- compteur + 1
        paramparallel[[c]][[i]] <- param[[compteur]]        
      } 
    }
    reference <- mclapply(X = paramparallel,
                          FUN = SpaTimeClusCpp,
                          input=input,
                          matT=matT, 
                          toollogistic = toollogistic,
                          mc.cores = nbcores, mc.preschedule = TRUE, mc.cleanup = TRUE)
    loglike <- rep(-Inf, nbcores)
    for (c in 1:nbcores) loglike[c] <- reference[[c]]@criteria@loglike
    reference <- reference[[which.max(loglike)]]
  }else{
    reference <- SpaTimeClusCpp(input, param,  matT, toollogistic)
  }
  return(TuneOutput(reference))
}

###################################################################################
##' This function performs the maximum likelihood estimation for a known model in clustering
##'
##' 
##' @param ech \linkS4class{STCdata} It contains the observations to cluster (mandatory). 
##' @param G numeric. It defines possible numbers of components.
##' @param K numeric. It defines possible numbers of regressions per components
##' @param Q numeric. It defines possible degrees of regressions.
##' @param crit character. It indicates the criterion used for the model selection ("AIC", "BIC" or "ICL", optional, default is "BIC").
##' @param tol numeric. The algorithm is stopped when the loglikelihood increases less than tol during two successive iterations (optional, default is 0.1).
##' @param param list of \linkS4class{STCparam}. It gives the initial values of the EM algorithm (optional).
##' @param nbcores numeric.  It defines the numerber of cores used by the alogrithm, only for Linux and Mac (optional, default is 1).
##' @param nbinitSmall numeric. It defines the number of random initializations (optional, default is 100).
##' @param nbinitKept numeric. It defines the number of chains estimated until convergence (optional, default is 10).
##' @param nbiterSmall numeric. It defines the number of iterations before keeping the nbinitKept best chains (optional, default is 10).
##' @param nbiterKept numeric. It defines the maximum number of iterations before to stop the algorith; (optional, default is 500).
##' 
##'  
##' @return Returns an instance of \linkS4class{STCresults}.
##' @examples
##' data(ech)
##' 
##' output <- spatimeclus(ech = ech, G = 2, K = 2, Q = 2, crit = "BIC", nbcores = 6)
##' 
##' summary(output)
##' 
##' print(output)
##' 
##' plot(output)
##' 
##' 
##' @export
##'
##'
spatimeclus <- function(ech, G, K, Q, crit="BIC", tol=0.1, param=NULL, nbcores=1, nbinitSmall=100, nbinitKept=10, nbiterSmall=10, nbiterKept=500){ 
  if (nbinitSmall<nbinitKept) nbinitKept <- nbinitSmall
  listmodels <- list()
  for (g in G){
    for (k in K){
      for (q in Q){
        listmodels[[length(listmodels)+1]] <- STCmodel(g, k, q)
      }
    }
  }
  results <- list()
  allcrit <- rep(-Inf, length(results))
  for (it in 1:length(listmodels)){
    #cat("model ", it, "\n")
    results[[it]] <- spatimeclusModelKnown(ech, listmodels[[it]], tol, param, nbcores, nbinitSmall, nbinitKept, nbiterSmall, nbiterKept)
    if (crit=="BIC"){
      allcrit[it] <- results[[it]]@criteria@BIC
    }else if (crit=="AIC"){
      allcrit[it] <- results[[it]]@criteria@AIC      
    }else if (crit=="ICL"){
      allcrit[it] <- results[[it]]@criteria@ICL      
    }
  }
  return(results[[which.max(allcrit)]])
}


###################################################################################
##' This function performs classifies a sample 
##'
##' 
##' @param ech \linkS4class{STCdata} It contains the observations to classify (mandatory). 
##' @param model \linkS4class{STCmodel}. It defines the model at hand (mandatory).
##' @param param list of \linkS4class{STCparam}. It gives the model parameters (mandatory).
##' 
##'  
##' @return Returns an instance of \linkS4class{STCresults}.
##' @export
##'
##'
spatimeclass <- function(ech, model, param){  
  matT <- t(exp(sweep(log(matrix(1:ech@TT, model@Q+1, ech@TT, byrow=TRUE)), 1, 0:model@Q, "*")))
  toollogistic <- cbind(ech@map, rep(1, ech@JJ))
  for (tt in 2:ech@TT)    toollogistic <- rbind(toollogistic, cbind(ech@map, rep(tt,  ech@JJ)))
  toollogistic <- cbind(toollogistic, rep(1, nrow(toollogistic)))
  # Proba post computation
  proba <- Probacond(ech, model, param, matT, toollogistic)
  output<- new("STCresults", model=model, data=ech, param=param, criteria=new("STCcriteria", loglike= sum(log(rowSums(exp(sweep(proba$logcondinter, 1, apply(proba$logcondinter,1, max), "-")))) + apply(proba$logcondinter,1, max))))
  return(TuneOutput(output))
}