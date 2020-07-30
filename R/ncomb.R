#' Finds number of combinations for DVI problem
#'
#' The number of victims and missing persons and their sex is given. The number
#' of possible moves, i,e., the number of ways the victims can be identified
#' with the missing persons is calculated.
#' @param nVfemales Integer. Number of female victims.
#' @param nMPfemales Integer. Number of missing persons, females.
#' @param nVmales Integer. Number of male victims.
#' @param nMPmales Integer. Number of missing persons,  males.
#'
#' @return Number of possible moves
#' @export
#' @importFrom arrangements
#' @examples
#' # There are 4 female victims and four female missing persons
#' ncomb(4, 4, 0, 0) # = 209
#'
#' \dontrun{
#' # This can be check by counting,all females below
#' ids.to = c("MP1", "MP2", "MP3", "MP4") # Female victims
#' moves = list(V1 = c("V1", ids.to), V2 = c("V2", ids.to),
#'              V3 = c("V3", ids.to), V7 = c("V7", ids.to))
#' foo = expand.grid.nodup(moves)
#' length(foo)
#' }
ncomb = function(nVfemales, nMPfemales, nVmales, nMPmales){
   fmoves = function(V,M){
    nmoves = 0
    for (k in 0:min(V,M)){
      n1 = arrangements::ncombinations(V, k)
      n2 = arrangements::ncombinations(M, k)
      nmoves = nmoves + n1*n2*factorial(k)
    }
    nmoves
  }
  if(nVfemales < 0| nMPfemales < 0 | nVmales < 0 | nMPmales < 0)
    stop("All parameters must be non-negative integers")
  nfe = fmoves(nVfemales, nMPfemales) 
  nma = fmoves(nVmales, nMPmales)
  nfe*nma
}

