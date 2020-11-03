#' Finds number of assignments for DVI problem
#'
#' The number of victims and missing persons and their sex is given. The number
#' of possible assignments, i,e., the number of ways the victims can be identified
#' with the missing persons,
#' is calculated.
#' @param nVfemales Integer. Number of female victims.
#' @param nMPfemales Integer. Number of missing persons, females.
#' @param nVmales Integer. Number of male victims.
#' @param nMPmales Integer. Number of missing persons,  males.
#'
#' @return Number of possible assignments
#' @export
#' @importFrom arrangements ncombinations
#' @examples
#' # Example 1. 
#' # With three male victims, and two MP-s, that may or may not 
#' # belong to different families, the number of apriori possible solutions, 
#' # the number of `assignments` is
#' m1 = ncomb(0,0,3,2) # 13 
#' # Alternatively, we can generate the assignments
#' m2 = expand.grid.nodup(list(V1 = c("V1", "MP1", "MP2"), V2 = c("V2","MP1", "MP2"), 
#'                             V3 = c("V3", "MP1", "MP2")))
#' # and check that the number of assignments coincides
#' stopifnot(m1 == length(m2)) 
#' 
#' 
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

