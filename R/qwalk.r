#' @title Quantum Walk
#' @author Pan Binghuang
#' @description Generate the modes, the probabilities that the walker being found at vertices.
#' An adjacency matrix is need for the process.
#' @usage qwdap.qwalk(edges, startindex, lens, scals, getfloat)
#'
#' @param edges your N*N adjacency matrix saved as list.
#' @param startindex the initial position of the quantum walker.
#' @param lens the number of records required in a round of sampling by a scaling factor.
#' Set the length of the series according to requirements.
#' @param scals the scaling factors used.
#' @param getfloat Whether to return floating point data.
#'
#' @return a object of class 'CTQW', the quantum walk results and some parameters.
#' @details 'qwdap.qwalk()' is used to generated modes for time series analysis, the result is a object
#' of class 'CTQW', the modes are saved in the object as a 3-dim array, and the parameters are also
#' store in the object.
#' The continuous time quantum walk is a continuous process, the modes are generated with a series
#' of times, the parameter 'scals' can be understood as the tolerance of the arithmetic time series.
#' Multiply tolerances can be passed in to obtain modes on different time scales through parameter
#' 'scals'.
#' The probability of the series with the probabilities that the walker being found at the vertices, 
#' and the length depends on parameter 'lens'.
#' The data generated by this function is not recorded from the initial state. 
#' The shortest distance between all vertices and the initial position of the quantum walker is
#' obtained by the Dijkstra algorithm.The probabilities corresponding to each vertex are
#' recorded starting from the vertex furthest in the shortest distance is not 0.
#' The function is single thread.
#' 
#' @useDynLib QWDAP
#' @importFrom Rcpp evalCpp
#' @export qwdap.qwalk
#'
#' @examples
#' edges <- matrix(c(0,1,0,0,0,0,0,
#'                   1,0,1,0,0,0,0,
#'                   0,1,0,1,0,0,0,
#'                   0,0,1,0,1,0,0,
#'                   0,0,0,1,0,1,0,
#'                   0,0,0,0,1,0,1,
#'                   0,0,0,0,0,1,0),
#'                 nrow = 7)
#' res.qwalk <- qwdap.qwalk(edges,1,100,scals=seq(from=0.01, by=0.01, length.out=5))
#' 
qwdap.qwalk<-function(edges, startindex = 1, lens=100, scals=c(0.01), getfloat=FALSE){
  if(!is.matrix(edges)){
    stop("The parameter 'edges' is not a matrix.")
  }
  if(nrow(edges)!=ncol(edges)){
    stop("The matrix 'edges' is not a square matrix.")
  }
  if(!is.vector(scals) || !is.numeric(scals)){
    stop("The parameter 'scals' is not a numeric vector.")
  }
  if(startindex<1 || startindex>nrow(edges)){
    stop("The parameter 'startindex' is out of range.")
  }
  if(!length(scals)){
    stop("The length of 'scals' is zero.")
  }
  
  dijkstra <- function(mgraph, v0){
    mgraph <- as.data.frame(mgraph)
    if(v0 < 1 || v0 > length(mgraph)){
      stop("v0 is an error value.")
    }
    for(i in c(1:nrow(mgraph))){
      for(j in c(1:ncol(mgraph))){
        if(i==j){
          mgraph[i,j] = 0
          next
        }
        if(mgraph[i,j] == 0){
          mgraph[i,j] = Inf
        }
      }
    }
    N = length(mgraph)
    final = p = d = as.data.frame(matrix(rep(0,length(mgraph)),nrow = 1))
    colnames(p) = colnames(d) = colnames(mgraph)
    v = w = k = min1 = NA
    for(v in c(1:N)){
      final[[v]] = 0
      d[[v]] = mgraph[v0,v]
      if(mgraph[v0,v] < Inf){
        p[[v]] = v0
      }else{
        p[[v]] = -1
      }
    }
    p[[v0]] = -1
    final[[v0]] = 1
    for(v in c(2:N)){
      if(v < 2) break
      min1 = Inf
      for(w in c(1:N)){
        if(!final[[w]] && d[[w]] < min1){
          k = w
          min1 = d[[w]]
        }
      }
      final[[k]] = 1
      for(w in c(1:N)){
        if(!final[[w]] && ((min1 + mgraph[k,w]) < d[[w]])){
          d[[w]] = min1 + mgraph[k,w]
          p[[w]] = k
        }
      }
    }
    return(list("dis" = d,"path" = p))
  }
  
  dij_edges = dijkstra(edges,startindex)
  Tag_node = which(dij_edges$dis==max(dij_edges$dis),arr.ind=TRUE)[1,2]
  
  N = nrow(edges)
  times = length(scals)
  
  ctqw = qwalkRcpp(edges,startindex-1,lens,scals,Tag_node-1,getfloat,100)
  
  # if(!getfloat){
  #   ctqw = round(ctqw*100)
  # }
  
  ar_ctqw = array(dim = c(lens, N, times))
  
  pointer = 0
  for(i in c(1:times)){
    ar_ctqw[,,i] = ctqw[(pointer+1):(pointer+lens),]
    pointer = pointer + lens
  }
  
  # if(is.null(alldata)){
  #   alldata <- array(t_data, dim = c(dim(t_data),i),
  #                    dimnames = list(c(1:nrow(t_data)),colnames(edges),paste("S",c(1:i),sep="")))
  # }else{
  #   alldata <- array(c(alldata, t_data), dim = c(dim(t_data),i),
  #                    dimnames = list(c(1:nrow(t_data)),colnames(edges),paste("S",c(1:i),sep="")))
  # }

  # res data.frame
  res<-list(ctqw = ar_ctqw,edges=edges,startindex=startindex,
            scals = scals, lens=lens, getfloat=getfloat, tag_node=Tag_node)
  res<-structure(res,class="CTQW")
  return(res)
}
