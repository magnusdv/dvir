#' LR for discrete property in DVI problems
#'
#' Based 
#'
#' @param a Character vector. An assignment
#' @param x List with 1 (has property) or 0 (do not have property) for victims
#' @param y List with 1 (has property) or 0 (do not have property) for missing persons
#' @param alpha Double. P(Y=1)
#' @param mu Double. P(X = 1 | Y = 0) = P(X = 0 | Y = 1)
#' @details The LR 
#'
#' @return LR
#' @seealso [expand.grid.nodup2()]
#' @examples
#' # Example 1
#' miss = c('*', example2$missing)
#' lst = list(V1 = miss, V2 = miss, V3 = miss)
#' tab = expand.grid.nodup2(lst)
#' nV = dim(tab)[2]
#' a = tab[1, 1:nV]
#' a[1:nV] = '*'
#' x = list(V1 = 1, V2 = 0, V3 = 0)
#' y = list(M1 = 1, M2 = 0, M3 = 0)
#' LRDiscrete(a, x, y, alpha = 0.1, mu = 0.01)
#' na = dim(tab)[1]
#' foo1 = rep(NA, na)
#' for (i in 1:na)
#' foo1[i] = LRDiscrete(tab[i,1:nV], x, y, alpha = 0.05, mu = 0.01)
#' res = data.frame(tab, foo1)
#' res = res[order(res$foo1, decreasing = T),]
#' head(res)
#' 
#' @export

LRDiscrete = function(a, x, y, alpha = 0.5, mu = 0.05){
  nV = length(a)
  index = (1:nV)[a != '*']
  n00 = n10 = n01 = n11 = 0
  victims = names(a)
  for (i in index){
   V = victims[i]
   x1 = as.double(x[[V]])
   M = a[[i]]
   y1 = as.double(y[[M]])
   if (x1 == 0 & y1 == 0)
     n00 = n00 + 1
   else if (x1 == 1 & y1 == 0)
     n10 = n10 + 1
   else if (x1 == 0 & y1 == 1)
     n01 = n01 + 1
   else
     n11 = n11 + 1
   }

   LR00 = (1-mu)/(1-alpha)
   LR01 = mu/(1-alpha)
   LR10 = mu/alpha
   LR11 = (1-mu)/alpha
   LR = LR00^n00 * LR01^n01 * LR10^n10 * LR11^n11
   LR
}
