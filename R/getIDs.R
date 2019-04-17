#' Get IDs from list of singeltons or pedigrees
#'
#' @param peds list of sigletons or pedigree
#' @export
#' @return Character vector with IDs
#' @examples 
#' pm1 = singleton("PM1")
#' getIDs(pm1)
#' pm2 = singleton("PM2")
#' getIDs(list(pm1, pm2))
#' fam1 = nuclearPed(1, father = "R1", mother = "R2", child= "MP1" )
#' getIDs(fam1)
#' fam2 = nuclearPed(1, father = "R3", mother = "R4", child= "MP2" )
#' am = list(fam1, fam2)
#' getIDs(am)
getIDs = function(peds){
  if(is.singleton(peds) | is.ped(peds))
    return(peds$ID)
  if(is.pedList(peds) | all(unlist(lapply(peds, function(x) is.ped(x)))))
    return(unlist(lapply(peds, function(peds) peds$ID)))
}