#' Generate list of possible moves
#'
#'A a sex consistent list of moves for `global` is generated
#' @aliases generateMoves
#' @param from A list of singletons.
#' @param to A list of pedigrees.
#' @param ids.to Character vector with names of missing persons.
#' @details Identity moves are included.
#'
#' @return A list of moves 
#' @seealso `global`.
#' @import pedtools
#' @examples
#' \dontrun{
#' from = list(singleton("V1", 1), singleton("V2", 2))
#' ids.to = paste("MP", 1:4, sep = "")
#' to = list(nuclearPed(children = ids.to[1:3]),
#'           nuclearPed(children = ids.to[4], sex = 2))
#' generateMoves(from, to,  ids.to)
#' }
#' @export generateMoves
generateMoves = function(from, to,  ids.to){
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
  ids.from = as.character(lapply(from, function(x) x$ID))
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

