initparam <- function(ech, model){
  lambda <- beta <- sigma <- list()
  prop <- runif(model@G)
  prop <-  prop / sum(prop)
  for (g in 1:model@G){
    lambda[[g]] <- matrix(0, model@K, 4)
    if (model@K>1) lambda[[g]][-1,] <- matrix(runif((model@K-1)*4,min = -1,max = 1), model@K-1, 4)*0
    colnames(lambda[[g]]) <- c("row", "col", "time", "cste")
    rownames(lambda[[g]]) <- paste("polynom", 1:model@K,sep=".")
    beta[[g]] <- matrix(runif(model@K*(model@Q+1),min = -1,max = 1), model@K, model@Q+1)
    if (model@Q>0){
      for (k in 1:model@K){
        borne <- sort(sample(1:ech@TT, 2))
        smalltmp <- cbind(ech@x[sample(1:ech@JJ, 1) + ((borne[1]:borne[2])-1)*ech@JJ ,sample(1:ech@n, 1)], borne[1]:borne[2])
        if (model@Q>1){for (q in 2:model@Q)    smalltmp <- cbind(smalltmp, (borne[1]:borne[2])**q)}
        beta[[g]][k,] <- c(mean(ech@x[sample(1:ech@JJ, 1) + ((borne[1]:borne[2])-1)*ech@JJ ,sample(1:ech@n, 1)]), rep(0, model@Q)) 
      }
    }
  }
  sigma <- matrix(max(var(ech@x[,1])), model@G, model@K)
  rownames(sigma) <- names(lambda) <- names(beta) <- names(prop) <- paste("compo", 1:model@G, sep=".")
  colnames(sigma) <- paste("polynom", 1:model@K, sep=".")
  return(new("STCparam", proportions=prop, lambda=lambda, beta=beta, sigma=sigma))
}

Probacond <- function(ech, model, param, matT, toollogistic){
  output <- list(condintra=list(), condintramargin=list(), logcondinter=matrix(0, ech@n, model@G, byrow = TRUE))
  for (g in 1:model@G){
    output$condintra[[g]] <-  list()
    output$condintramargin[[g]] <- matrix(0, ech@JJ*ech@TT, ech@n)
    poidspolynom <- matrix(0, ech@TT*ech@JJ, model@K)
    for (k in 1:model@K) poidspolynom[,k] <- toollogistic %*% param@lambda[[g]][k,]
    poidspolynom<- exp(sweep(poidspolynom,1,apply(poidspolynom,1,max), "-"))
    poidspolynom <- poidspolynom / rowSums(poidspolynom)
    
    for (k in 1:model@K){
      output$condintra[[g]][[k]] <- matrix(poidspolynom[,k], ech@JJ*ech@TT, ech@n) * dnorm( sweep(ech@x, 1, rep(as.matrix(matT[, 1:(model@Q+1)]) %*% (param@beta[[g]][k,]), each=ech@JJ), "-"), sd = sqrt(param@sigma[g,k]))
      output$condintramargin[[g]] <- output$condintramargin[[g]] + output$condintra[[g]][[k]]
    }  
    output$logcondinter[,g] <- colSums(log(output$condintramargin[[g]])) + log(param@proportions[g]) 
  }
  return(output)
}


TuneOutput <- function(output){
  output@param@proportions <- as.numeric(output@param@proportions)
  for (g in 1:output@model@G){
    output@param@lambda[[g]] <- matrix(output@param@lambda[[g]], output@model@K, 4)
    colnames(output@param@lambda[[g]]) <- c("coordinate1", "coordinate2", "time", "intercept")
    output@param@beta[[g]] <- matrix(output@param@beta[[g]], output@model@K, output@model@Q+1)
    colnames(output@param@beta[[g]]) <- paste("degree", 0:output@model@Q, sep=".")
    rownames(output@param@lambda[[g]]) <- rownames(output@param@beta[[g]]) <- paste("polynom", 1:output@model@K, sep=".")
  }
  names(output@param@beta) <- names(output@param@lambda) <- paste("component", 1:output@model@G, sep=".")
  
  matT <- t(exp(sweep(log(matrix(1:output@data@TT, output@model@Q+1, output@data@TT, byrow=TRUE)), 1, 0:output@model@Q, "*")))
  toollogistic <- cbind(output@data@map, rep(1, output@data@JJ))
  for (tt in 2:output@data@TT)    toollogistic <- rbind(toollogistic, cbind(output@data@map, rep(tt,  output@data@JJ)))
  toollogistic <- cbind(toollogistic, rep(1, nrow(toollogistic)))
  proba <- Probacond(output@data, output@model, output@param, matT, toollogistic)
  sig <- weightlogistic <-list()
  for (g in 1:output@model@G){
    sig[[g]] <- list()
    weightlogistic[[g]]  <- matrix(0, output@data@JJ*output@data@TT, output@model@K)
    for (k in 1:output@model@K) names(output@param@beta[[g]][k,]) <- paste("coeff", 0:(length(output@param@beta[[g]][k,])-1), sep=".")
  }  
  
  output@partitions@hardseg <- list()
  for (g in 1:output@model@G){
    output@partitions@hardseg[[g]] <- matrix(apply(exp(toollogistic %*% t(output@param@lambda[[g]])), 1, which.max), output@data@JJ, output@data@TT)
    colnames(output@partitions@hardseg[[g]]) <- paste("time", 1:output@data@TT, sep=".")
    rownames(output@partitions@hardseg[[g]]) <- paste("site", 1:output@data@JJ, sep=".")
  }
  names(output@partitions@hardseg) <- paste("comp", 1:output@model@G, sep=".")
  
  
  output@partitions@fuzzyind <- exp(sweep(proba$logcondinter, 1, apply(proba$logcondinter,1, max), "-"))
  output@partitions@fuzzyind  <- output@partitions@fuzzyind  / rowSums(output@partitions@fuzzyind )
  output@partitions@hardind <- apply(output@partitions@fuzzyind, 1, which.max) 
  names(output@partitions@hardind) <- rownames(output@partitions@fuzzyind) <- colnames(output@data@x)
  colnames(output@partitions@fuzzyind) <- paste("comp", 1:output@model@G, sep=".")
  
  output@criteria@AIC <- output@criteria@loglike - output@model@nbparam
  output@criteria@BIC <- output@criteria@loglike - 0.5 * output@model@nbparam * log(output@data@n)
  rownames(output@param@sigma) <- paste("component", 1:output@model@G, sep=".")
  colnames(output@param@sigma) <- paste("polynom", 1:output@model@K, sep=".")
  entrop <- 0
  for (g in 1:output@model@G) entrop <- entrop + sum(log(output@partitions@fuzzyind[which(output@partitions@hardind==g),g]))
  output@criteria@ICL <- output@criteria@BIC + entrop
  return(output)
}