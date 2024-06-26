#' @title RReliefF
#' @description Mode selection by RReliefF. The purpose of this function
#' is to select the part modes with similar characteristics to the observed time series
#' from the modes generated by the quantum walk. And it is based on the data model.
#' @usage qwdap.rrelieff(real, ctqw, index, num, plotting)
#' @param real the real series observed.
#' @param ctqw the 'CTQW' object.
#' @param index the index of the data for mode selection.
#' @param num the number of series required.
#' @param plotting whether to plot.
#'
#' @return a 'QWMS' object.
#' 
#' @details The 'QWMS' object include the original time series and the modes generated by
#' quantum walks.
#' @import CORElearn
#' @export qwdap.rrelieff
#'
#' @examples
#' data("traffic.qw")
#' data("trafficflow")
#' res.rrelieff <- qwdap.rrelieff(trafficflow,traffic.qw,1,30,TRUE)
#' 
qwdap.rrelieff<-function(real,ctqw,index,num = -1,plotting = FALSE){
  # library(CORElearn)
  # get the data
  if(!inherits(ctqw, 'CTQW')){
    stop("The parameter 'ctqw' is not a 'CTQW' object.")
  }
  if(nrow(real)!=ctqw$lens){
    stop("The row of 'real' is not equal to the 'lens' of 'ctqw'.")
  }
  if(!is.data.frame(real)){
    real = as.data.frame(real)
  }
  proc.data <- cbind(real[index],as.data.frame(ctqw$ctqw[,index,]))
  res <- attrEval(1,proc.data,estimator = "RReliefFbestK")
  res <- res[order(abs(res),decreasing = T)]
  res <- t(data.frame(res))
  rownames(res)<-"importance"
  
  if(plotting){
    plot(x=c(1:ncol(res)),y=abs(res[1,]),type="p",lwd=2,xlab="",ylab="")
  }
  if(num != -1){
    res <- subset(res, select = colnames(res)[1:num])
  }
  res<-list(real=as.matrix(real[index]), ctqw=ctqw$ctqw[,index,], index = index,
            method = "RReliefFbestK", variate = colnames(res),importance=res)
  res<-structure(res,class="QWMS")
  return(res)
}
