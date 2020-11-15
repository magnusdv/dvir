#' The number of assignments for DVI problem
#'
#' The number of victims and missing persons of each sex is given. The number
#' of possible assignments, i,e., the number of ways the victims can be
#' identified with the missing persons, is calculated.
#'
#' @param nVfemales Integer. The number of female victims.
#' @param nMPfemales Integer. The number of female missing persons.
#' @param nVmales Integer. The number of male victims.
#' @param nMPmales Integer. The number of male missing persons.
#'
#' @return The total number of possible assignments.
#'
#' @examples
#'
#' # Example: 3 male victims; 2 male missing persons.
#' # The number of a priori possible assignments is
#' m1 = ncomb(0,0,3,2) # 13
#'
#' # The complete list of assignments
#' m2 = expand.grid.nodup(list(V1 = c("V1", "MP1", "MP2"),
#'                             V2 = c("V2", "MP1", "MP2"),
#'                             V3 = c("V3", "MP1", "MP2")))
#'
#' # Check that the number of assignments coincides
#' stopifnot(m1 == length(m2))
#'
#' @export
ncomb = function(nVfemales, nMPfemales, nVmales, nMPmales){
  fmoves = function(V,M){
    nmoves = 0
    for (k in 0:min(V,M)){
      n1 = choose(V, k)
      n2 = choose(M, k)
      nmoves = nmoves + n1*n2*factorial(k)
    }
    nmoves
  }
  
  if(nVfemales < 0| nMPfemales < 0 | nVmales < 0 | nMPmales < 0)
    stop("All parameters must be non-negative integers")
  
  nfe = fmoves(nVfemales, nMPfemales) 
  nma = fmoves(nVmales, nMPmales)
  nfe * nma
}

