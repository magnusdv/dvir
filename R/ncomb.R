#' Finds number of combinations for DVI problem
#' 
#' A formula is implemented
#' @param Vfe Integer. Number of female victims
#' @param Mma Integer. Number of male victims
#' @param Vfe Integer. Number of female missing persons
#' @param Vma Integer. Number of male missing persons
#' 
#' @return Number of possible moves
#' @export
#' @examples 
#' There 3 female victims, 2 male victims,6 missing females, 6 missing males
#' #' library(arrangements)
#' ncomb(3,6,2,6) # = 9847 As checked by counting
#' for nfi example:
#' \dontrun{
#' moves =list(c(1,3,8,12,5,9,101),
#' c(7,2,10,11,6,4,102),
#' c(3,12,1,8,5,9,103),
#' c(4,7,2,10,11,6,104),
#' c(5,12,3,8,9,1,105))
#' foo = expand.grid.nodup(moves)
#' length(foo)
#' }
ncomb = function(Vfe, Mfe, Vma, Mfa){
  fmoves = function(V,M){
    nmoves = 0
    for (k in 0:min(V,M)){
      n1 = ncombinations(V, k)
      n2 = ncombinations(M, k)
      nmoves = nmoves + n1*n2*factorial(k)
    }
    nmoves
  }
  if(Vfe < 0| Mfe < 0 | Vma < 0 | Mfa < 0)
    stop("All parameters must be non-negative integers")
  nfe = fmoves(Vfe, Mfe) 
  nma = fmoves(Vma, Mfa)
  nfe*nma
}

