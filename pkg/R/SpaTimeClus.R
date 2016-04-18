##' SpaTimeClus a package for clustering spatio-temporal data
##'
##' SpaTimeClus is a tool for clustering Spatio-Temporal data.
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
##' Author: Cheam A., Marbac M., and McNicholas P.
##'
##' @references Cheam A., Marbac M., and McNicholas P., Model-Based Clustering for Spatio-Temporal Data Applied for Air Quality.
##' 
##' @examples
##' data(airparif)
##' 
##' # Clustering of the data by considering the spatial dependencies
##' res.spa <- spatimeclus(airparif$obs,  G=3, K=4, Q=4, map = airparif$map, nbinitSmall=50, nbinitKept=5, nbiterSmall=5)
##' summary(res.spa)
##' 
##' # Clustering of the data without considering the spatial dependencies
##' res.nospa <- spatimeclus(airparif$obs,  G=3, K=4, Q=4, nbinitSmall=50, nbinitKept=5, nbiterSmall=5)
##' summary(res.nospa)
##' 
NULL



spatimeclusModelKnown <- function(obs, model, param=NULL, tune=tune){
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
    for (it in 1:tune@nbinitSmall)  param[[it]] <- initparam(obs, model, matT)
  }
  input <-new("STCresults", 
              model=model, 
              data=obs, 
              criteria=new("STCcriteria", loglike= -Inf),
              tune=tune
  )  
  return(TuneOutput(SpaTimeClusCpp(input, param,  matT) ))
}

###################################################################################
##' This function performs the maximum likelihood estimation for a known model in clustering
##'
##' @param x array It contains the observations to cluster where the dimesions are respectively: number of the observation, site of the observation, time of the observation.  
##' @param G numeric. It defines possible numbers of components.
##' @param K numeric. It defines possible numbers of regressions per components
##' @param Q numeric. It defines possible degrees of regressions.
##' @param map matrix. It gives the spatial coordiantes of each site.
##' @param m numeric. It indicates the moments of observations.
##' @param crit character. It indicates the criterion used for the model selection ("AIC", "BIC" or "ICL", optional, default is "BIC").
##' @param tol numeric. The algorithm is stopped when the loglikelihood increases less than tol during two successive iterations (optional, default is 0.01).
##' @param param list of \linkS4class{STCparam}. It gives the initial values of the EM algorithm (optional).
##' @param nbcores numeric.  It defines the numerber of cores used by the alogrithm, only for Linux and Mac (optional, default is 1).
##' @param nbinitSmall numeric. It defines the number of random initializations (optional, default is 500).
##' @param nbinitKept numeric. It defines the number of chains estimated until convergence (optional, default is 30).
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
spatimeclus <- function(obs, G, K, Q, map=NULL, m=1:(dim(obs)[3]), crit="BIC", tol=0.001, param=NULL, nbcores=1, nbinitSmall=500, nbinitKept=50, nbiterSmall=20, nbiterKept=500){
  obs <- BuildSTCdata(obs, map, m=m)
  nbcores <- min(detectCores(all.tests = FALSE, logical = FALSE),  nbcores)
  if (nbinitSmall<nbinitKept) nbinitKept <- nbinitSmall
  listmodels <- list()
  resallmodels <- matrix(NA, length(G)*length(K)*length(Q),3)
  lig <- 0
  for (g in G){
    for (k in K){
      for (q in Q){
        listmodels[[length(listmodels)+1]] <- STCmodel(g, k, q, is.null(map))
        lig <- lig + 1
        resallmodels[lig,] <- c(g, k, q)
      }
    }
  }
  for (i in 1:obs@n) obs@x[i,,-1] <- obs@x[i,,-1] - obs@x[i,,-obs@TT]
  tune <- new("STCtune", tol=tol, nbinitSmall=nbinitSmall, nbinitKept=nbinitKept, nbiterSmall=nbiterSmall, nbiterKept=nbiterKept)
  if (nbcores>1){
    results <- mclapply(X = listmodels,
                        FUN = spatimeclusModelKnown,
                        obs=obs,
                        param=param,
                        tune=tune,                      
                        mc.cores = nbcores, mc.preschedule = TRUE, mc.cleanup = TRUE)
  }else{
    results <- list()
    for (it in 1:length(listmodels))  results[[it]] <- spatimeclusModelKnown(obs, listmodels[[it]], param, tune)
  }
  allcrit <- rep(-Inf, length(results))  
  degen <- rep(-Inf, length(results))
  garde <- degen
  for (it in 1:length(listmodels)){
    degen[it] <- results[[it]]@criteria@degeneracy
    if (crit=="BIC"){
      allcrit[it] <- results[[it]]@criteria@BIC
    }else if (crit=="AIC"){
      allcrit[it] <- results[[it]]@criteria@AIC      
    }else if (crit=="ICL"){
      allcrit[it] <- results[[it]]@criteria@ICL      
    }
  }
  garde <- allcrit
  garde[which(degen>0.5)] <- -Inf
  results <- results[[which.max(garde)]]
  results@allmodels <- cbind(resallmodels, allcrit, degen)
  return(results)
}
