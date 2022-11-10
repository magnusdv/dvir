#' Sex-consistent pairings
#'
#' Generate a list of sex-consistent pairings for each victim in a DVI problem.
#' By default, the empty pairing (denoted `*`) is included for each victim.
#'
#' @param dvi A `dviData` object, typically created with [dviData()].
#' @param includeEmpty A logical. If TRUE (default), the do-nothing symbol (`*`)
#'   is included for each victim.
#' @param ignoreSex A logical.
#'
#' @return A list of character vectors. Each vector is a subset of `missing`,
#'   plus the character `*` denoting no pairing.
#'
#' @seealso [jointDVI()]
#'
#' @examples
#'
#' pm = list(singleton("V1", sex = 1),
#'           singleton("V2", sex = 2))
#'
#' missing = paste0("M", 1:4)
#' am = list(nuclearPed(children = missing[1:3]),
#'           nuclearPed(children = missing[4], sex = 2))
#' 
#' dvi = dviData(pm, am, missing)
#' generatePairings(dvi)
#'
#' @export
generatePairings = function(dvi, includeEmpty = TRUE, ignoreSex = FALSE){
  
  if(!inherits(dvi, "dviData"))
    stop2("First argument must be `dviData` object. (As of dvir version 3.0.0)")
  pm = dvi$pm
  am = dvi$am
  missing = dvi$missing
  
  # ID of victims
  vics = unlist(labels(pm), use.names = FALSE)
  
  if(ignoreSex) {
    mp = if(includeEmpty) c("*", missing) else missing
    lst = replicate(length(vics), mp, simplify = FALSE)
    names(lst) = vics
    return(lst)
  }
    
  # Vectors with sex of victims and missing
  sexVics = getSex(pm, vics, named = TRUE)
  sexMP = getSex(am, missing)
  
  # For each victim, find missing of same sex
  lst = lapply(vics, function(v) {
    mp = missing[sexMP == sexVics[v]]
    if(includeEmpty) 
      mp = c("*", mp)
    mp
  })
  names(lst) = vics

  lst
}

