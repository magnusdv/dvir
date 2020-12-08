#' Generate list of possible moves
#'
#' A sex-consistent list of moves for [global()] is generated
#'
#' @param from A list of singletons.
#' @param to A list of pedigrees.
#' @param ids.to Character vector with names of missing persons.
#' @param ignore.sex A logical.
#' @param expand.grid A logical. If TRUE, [expand.grid.nodup()] is applied to
#'   the result before returning.
#' @details Identity moves are included.
#'
#' @return A list of moves.
#' @seealso [global()]
#'
#' @examples
#'
#' library(pedtools)
#'
#' from = list(singleton("V1", sex = 1),
#'             singleton("V2", sex = 2))
#' ids.to = paste0("MP", 1:4)
#' to = list(nuclearPed(children = ids.to[1:3]),
#'           nuclearPed(children = ids.to[4], sex = 2))
#' generateMoves(from, to, ids.to)
#'
#' @import pedtools
#' @export
generateMoves = function(from, to, ids.to, ignore.sex = FALSE, expand.grid = FALSE){
  
  # ID of victims
  vics = unlist(labels(from), use.names = FALSE)
  
  if(ignore.sex) {
    lst = sapply(vics, function(v) c("*", ids.to), 
                 simplify = FALSE, USE.NAMES = TRUE)
    return(lst)
  }
    
  # Vectors with sex of victims amd MPs
  sexVics = getSex(from, vics, named = TRUE)
  sexMPs = getSex(to, ids.to)
  
  # For each victim, find MPs of same sex
  res = sapply(vics, function(v) c("*", ids.to[sexMPs == sexVics[v]]),
               simplify = FALSE, USE.NAMES = TRUE)

  if(expand.grid)
    res = expand.grid.nodup(res)
  
  res
}



generateMoves_old = function(from, to,  ids.to){
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
  ids.from = unlist(labels(from))
  nV = length(ids.from)
  sexVictims = getSex(from, ids.from)
  femaleVictims = ids.from[sexVictims == 2]
  maleVictims = ids.from[sexVictims == 1]
  nMP = length(ids.to)
  sexMPs = getSex(to, ids.to)
  femaleMPs = ids.to[sexMPs == 2]
  maleMPs =ids.to[sexMPs == 1]
  lstFemales = generate1(femaleVictims, femaleMPs)
  lstMales = generate1(maleVictims, maleMPs)
  c(lstFemales, lstMales)
}

