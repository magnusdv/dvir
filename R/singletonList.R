#' Creates a list of singlelots
#' 
#' Wraper for singleton function of pedtools
#' 
#' @param n integer Noof singletons
#' @param ids Character vector. Names of victims
#' @param famid Character vector. Names of families
#' @param sex Integer vector
#' @return  Alist of singletons
#' @examples 
#' singletonList()
#' @export
singletonList = function(n = 2, ids = paste("V", 1:n, sep = ""), 
                         famid = paste("F", 1:n, sep = ""), sex = rep(1,n)){
  singles = list()
  for (i in 1:n)
    singles[[i]] = singleton(ids[i], sex = sex[i], famid = famid[i])
  singles
}
