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
##' data(obs)
##' 
NULL



spatimeclusModelKnown <- function(obs, model, tol=0.1, param=NULL, nbcores=1, nbinitSmall=100, nbinitKept=10, nbiterSmall=10, nbiterKept=500){
  nbcores <- min(detectCores(all.tests = FALSE, logical = FALSE),  nbcores)
  nbinitSmall <- nbcores * ceiling(nbinitSmall / nbcores) 
  nbinitKept <- nbcores * ceiling(nbinitKept / nbcores) 

  if (model@Q>0){
    matT <- matrix(obs@m, obs@TT, 1)
    if (model@Q>1) {for (q in 2:model@Q) matT <- cbind(matT, obs@m ** q)}
    matT <- cbind(rep(1, obs@TT), matT)
  }else{
    matT <- matrix(1, obs@TT, 1)
  }
  matT[-1,] <- matT[-1,] - matT[-obs@TT,] 
  if (is.null(param)){
    param <- list()
    for (it in 1:nbinitSmall)  param[[it]] <- initparam(obs, model, matT)
  }
  input <-new("STCresults", 
              model=model, 
              data=obs, 
              criteria=new("STCcriteria", loglike= -Inf),
              tune=new("STCtune", 
                       tol=tol, 
                       nbinitSmall=ceiling(nbinitSmall / nbcores), 
                       nbinitKept=ceiling(nbinitKept / nbcores), 
                       nbiterSmall=nbiterSmall, 
                       nbiterKept=nbiterKept
                       )
              )  
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
                          matT=newtool, 
                          mc.cores = nbcores, mc.preschedule = TRUE, mc.cleanup = TRUE)

    loglike <- rep(-Inf, nbcores)
    for (c in 1:nbcores) loglike[c] <- reference[[c]]@criteria@loglike
    reference <- reference[[which.max(loglike)]]
  }else{
    print("deb")
    
    reference <- SpaTimeClusCpp(input, param,  matT)
    print("fin")
    
  }
  #return(TuneOutput(reference))
  return(reference)
}

###################################################################################
##' This function performs the maximum likelihood estimation for a known model in clustering
##'
##' 
##' @param obs \linkS4class{STCdata} It contains the observations to cluster (mandatory). 
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
##' data(obs)
##' 
##' 
##' 
##' @export
##'
##'
spatimeclus <- function(obs, G, K, Q, crit="BIC", tol=0.1, param=NULL, nbcores=1, nbinitSmall=100, nbinitKept=10, nbiterSmall=10, nbiterKept=500){ 
  if (nbinitSmall<nbinitKept) nbinitKept <- nbinitSmall
  listmodels <- list()
  for (g in G){
    for (k in K){
      for (q in Q){
        listmodels[[length(listmodels)+1]] <- STCmodel(g, k, q)
      }
    }
  }
  for (i in 1:obs@n) obs@x[i,,-1] <- obs@x[i,,-1] - obs@x[i,,-obs@TT]
  results <- list()
  allcrit <- rep(-Inf, length(results))
  for (it in 1:length(listmodels)){
    #cat("model ", it, "\n")
    results[[it]] <- spatimeclusModelKnown(obs, listmodels[[it]], tol, param, nbcores, nbinitSmall, nbinitKept, nbiterSmall, nbiterKept)
    if (crit=="BIC"){
      allcrit[it] <- results[[it]]@criteria@BIC
    }else if (crit=="AIC"){
      allcrit[it] <- results[[it]]@criteria@AIC      
    }else if (crit=="ICL"){
      allcrit[it] <- results[[it]]@criteria@ICL      
    }
  }
  #return(results[[which.max(allcrit)]])
  return(results[[1]])
}


###################################################################################
##' This function performs classifies a sample 
##'
##' 
##' @param obs \linkS4class{STCdata} It contains the observations to classify (mandatory). 
##' @param model \linkS4class{STCmodel}. It defines the model at hand (mandatory).
##' @param param list of \linkS4class{STCparam}. It gives the model parameters (mandatory).
##' 
##'  
##' @return Returns an instance of \linkS4class{STCresults}.
##' @export
##'
##'
spatimeclass <- function(obs, model, param){  

  
  toollogistic <- cbind(obs@map, rep(1, obs@JJ))
  for (tt in 2:obs@TT)    toollogistic <- rbind(toollogistic, cbind(obs@map, rep(tt,  obs@JJ)))
  if (model@Q>1){
    newtool <- toollogistic
    if (model@Q>1){for (q in 2:model@Q) newtool <- cbind(newtool,newtool[,1:3]**q)}
    newtool <- cbind(rep(1, nrow(newtool)), newtool)
    colnames(newtool) <- c("cste",paste(rep(c("spa1","spa2","tps"),model@Q), ".deg", rep(1:model@Q, each=3), sep=""))  
  }else{
    newtool <- matrix(1, nrow(toollogistic), 1)
  }
  toollogistic <- cbind(toollogistic, rep(1, nrow(toollogistic)))
  # Proba post computation
  proba <- Probacond(obs, model, param, newtool, toollogistic)
  output<- new("STCresults", model=model, data=obs, param=param, criteria=new("STCcriteria", loglike= sum(log(rowSums(exp(sweep(proba$logcondinter, 1, apply(proba$logcondinter,1, max), "-")))) + apply(proba$logcondinter,1, max))))
  return(TuneOutput(output))
}