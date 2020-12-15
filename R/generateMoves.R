#' Generate list of possible moves
#'
#' A sex-consistent list of moves for [jointDVI()] is generated
#'
#' @param pm A list of singletons.
#' @param am A list of pedigrees.
#' @param missing Character vector with names of missing persons.
#' @param ignore.sex A logical.
#' @param expand.grid A logical. If TRUE, [expand.grid.nodup()] is applied to
#'   the result before returning.
#' @details Identity moves are included.
#'
#' @return A list of moves.
#' @seealso [jointDVI()]
#'
#' @examples
#'
#' library(pedtools)
#'
#' pm = list(singleton("V1", sex = 1),
#'           singleton("V2", sex = 2))
#' 
#' missing = paste0("M", 1:4)
#' am = list(nuclearPed(children = missing[1:3]),
#'           nuclearPed(children = missing[4], sex = 2))
#' generateMoves(pm, am, missing)
#'
#' @import pedtools
#' @export
generateMoves = function(pm, am, missing, ignore.sex = FALSE, expand.grid = FALSE){
  
  # ID of victims
  vics = unlist(labels(pm), use.names = FALSE)
  
  if(ignore.sex) {
    lst = sapply(vics, function(v) c("*", missing), 
                 simplify = FALSE, USE.NAMES = TRUE)
    return(lst)
  }
    
  # Vectors with sex of victims amd missing
  sexVics = getSex(pm, vics, named = TRUE)
  sexMP = getSex(am, missing)
  
  # For each victim, find missing of same sex
  res = sapply(vics, function(v) c("*", missing[sexMP == sexVics[v]]),
               simplify = FALSE, USE.NAMES = TRUE)

  if(expand.grid)
    res = expand.grid.nodup(res)
  
  res
}



generateMoves_old = function(pm, am, missing){
  # Function to generate for one sex
  generate1 = function(victims, missing){
    nV = length(victims)
    nMP = length(missing)
    if(nV < 1)
      return(NULL)
    lst = rep(list(missing), length(victims))
  
    # Add in identity moves
    for (i in 1:nV)
      lst[[i]] = c(victims[i], lst[[i]])
    names(lst) = victims
    lst
  }
  
  # Find female and male victims and MP-s
  ids.pm = unlist(labels(pm))
  nV = length(ids.pm)
  sexVictims = getSex(pm, ids.pm)
  femaleVictims = ids.pm[sexVictims == 2]
  maleVictims = ids.pm[sexVictims == 1]
  nMP = length(missing)
  sexMP = getSex(am, missing)
  femaleMP = missing[sexMP == 2]
  maleMP = missing[sexMP == 1]
  lstFemales = generate1(femaleVictims, femaleMP)
  lstMales = generate1(maleVictims, maleMP)
  c(lstFemales, lstMales)
}

