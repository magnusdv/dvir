#' Find single candidate list
#' 
#' @param vp Character vector with names of victims
#' @param pm Character vector with names of missing persons
#' @param sex list of sex, one component for victims, 
#' one for missing persons
#' @return Matrix of candidate moves
#' @examples 
#' vp = c("V1", "V2")
#' mp = c("MP1", "MP2")
#' mat1 = rbind(c("V1", 1),c("V2", 2))
#' dimnames(mat1) = list(vp, c("id", "sex"))
#' mat2 = rbind(c("MP1", 1),c("MP2", 2))
#' dimnames(mat2) = list(mp, c("id", "sex"))
#' sex = list(pmsex = mat1, amsex = mat2)
#' loop2(c("V1", "V2"), c("MP1", "MP2"), sex)
#' 
#' 
#' @export
loop2 <-
function(vp, mp, sex){
    loop1 = loop(vp, mp)
    d1 = dim(loop1)[1]
    index = rep(NA, d1)
    for (i in 1:d1){
      from1 = as.character(loop1[i,1])
      to1 = as.character(loop1[i,2])
      index[i] = sex[[1]][from1, 2] == sex[[2]][to1, 2]
    }
    ni = sum(index)
    if (ni == 0) return(NULL)
    loop1 = matrix(loop1[index, ], ncol = 2)
    colnames(loop1) = c("from", "to")
    loop1
  }

loop <-
function(from, to){
    from2 = rep(from, each = length(to))
    to2 = rep(to, length(from))
    loop1 = cbind(from2,to2)
    loop1
  }
