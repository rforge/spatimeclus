#'
#' plot function.
#' 
#' This function gives the plot of an instance of \code{\linkS4class{STCresults}}.
#' 
#' @param object instance of \code{\linkS4class{STCresults}}.
#' 
#' @name plot
#' @rdname plot-methods
#' @docType methods
#' @exportMethod plot
#' 
#' 
NULL

#' @rdname plot-methods
#' @aliases plot plot,STCresults-method
#'
#'
setMethod(
  f="plot",
  signature = c("STCresults"),
  definition = function(x,...){
    par(mfrow=c(x@model@G, 2))
    for (g in 1:x@model@G){
      plot(NA,  xlab="time", ylab="site", xlim=c(-0.3, 1), ylim=c(0,1), axes=FALSE)
      text(rep(-0.15, x@data@JJ), seq(0,1, length.out = x@data@JJ), rownames(x@data@map), cex = 0.5)
      image(t(x@partitions@hardseg[[g]]), col=1:x@model@K,add = TRUE)
      values <- matrix(0, x@model@K, x@data@TT)
      for (k in 1:x@model@K){
        mat <- rep(1, x@data@TT)
        if (x@model@Q!=0){
          for (q in 1:x@model@Q)        mat <- cbind(mat, (1:x@data@TT)**q)
        }
        values[k,] <- as.matrix(mat)%*%as.numeric(x@param@beta[[g]][k,])
      }
      plot(NA, xlab="time", ylab="value", xlim=c(1,x@data@TT), ylim=range(x@data@x[,which(x@partitions@hardind==g)]))
      for(i in which(x@partitions@hardind==g)){
        tmp <- matrix(x@data@x[,i] , x@data@JJ, x@data@TT)
        for (j in 1:x@data@JJ) lines(1:x@data@TT, tmp[j,], col="grey")
        
      }
      for (j in 1:x@data@JJ){
        tmp <- rep(0, x@data@TT)
        for (k in  1:x@model@K){
          tmp[which(x@partitions@hardseg[[g]][j,]==k)] <- values[k,which(x@partitions@hardseg[[g]][j,]==k)]
        }
        lines(1:x@data@TT, tmp)
      }
       
    }
  }
)
