#' @title Prediction
#' @description Based on the established model, make predict.
#' @usage qwdap.predict(in_model, data_range)
#'
#' @param in_model a 'QWMODEL' object, which is the model built by Stepwise Regression,
#' PCR, PLSR, PPR, VAR in this package.
#' @param data_range indicate the index range of the part data generated by quantum walks for predict.
#'
#' @return the predict data.
#' @import MTS
#' @export
#'
#' @examples
#' data(traffic.model.n1)
#' res.predict <- qwdap.predict(traffic.model.n1,c(501,720))
qwdap.predict <- function(in_model, data_range){
  if(class(in_model)!="QWMODEL"){
    print("The 'in_model' is not a 'QWMODEL' object.")
    return()
  }
  if(!is.vector(data_range)||!is.numeric(data_range)||length(data_range)<2||data_range[1]>data_range[2]){
    print("The 'data_range' is error.")
    return()
  }
  res<-NULL
  if(in_model$method == 'VAR'){
    res<-VARXpred(in_model$model, newxt = in_model$ctqw[data_range[1]:data_range[2],],
                  hstep = data_range[2]-data_range[1]+1)
  }else if(in_model$method %in% c("Stepwise Regression", "PCR", "PLSR")){
    res<-predict(in_model$model, newdata = as.data.frame(in_model$ctqw[data_range[1]:data_range[2],]))
  }else if(in_model$method == "PPR"){
    res<-predict(in_model$model)
  }
  return(res)
}