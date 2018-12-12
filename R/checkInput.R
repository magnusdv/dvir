#checkInput
#' Checks input
#' 
#' @param pm A list of singeltons 
#' @param am A list of pedigrees
#' @param vp Character vector with names of victims
#' @param mp Character vector with names of missing persons
#' @return list with NULL if OK otherwise error message and null likelihood
#' @export
#' @examples
#' vp = c("V1", "V2")
#' pm = singletonList(2, ids = vp)
#' m = marker(pm[[1]], alleles = 1:2, "V1" = 1:2)
#' pm[[1]] = addMarkers(pm[[1]], m)
#' m = marker(pm[[2]], alleles = 1:2, "V2" = 1:2)
#' pm[[2]] = addMarkers(pm[[2]], m)
#' mp = c("MP1", "MP2")
#' am = nuclearPed(2, children = mp, father = "R1", mother = "R2")
#' m = marker(am, "R1" = 1:2, "R2"= 1:2)
#' am = addMarkers(am, m)
#' checkInput(pm, am, vp, mp)
#' checkInput(pm, am, vp[c(1,1)], mp)
#' checkInput(pm, am, vp, mp[c(1,1)])
#' m = marker(am, "R1" = 1:2, "R2"= 1:2, "MP1" = 1)
#' am = addMarkers(am, m)
#' checkInput(pm, am, vp, mp) 
#' am = setAlleles(am, alleles = 0)
#' m = marker(am, "R1" = 1:2, "R2"= 1:2)
#' am = addMarkers(am, m)
#' pm = setAlleles(pm[[1]], alleles = 0)
#' checkInput(pm, am, vp, mp)

checkInput = function(pm, am, vp , mp ){
  errortext = NULL
  if(!pedtools::is.pedList(pm) & !pedtools::is.ped(pm))
    errortext = c("First argument must be a  ped object or a list of such")
  if(!pedtools::is.pedList(am) & !pedtools::is.ped(am))
    errortext = c("Second argument must be a  ped object or a list of such")
  nvp = length(vp)
  if(nvp < 1) errortext = c("Third argument must be a vector of length 1 or more")
  nmp = length(mp)
  if(nmp < 1) errortext = c("Four argument must be a vector of length 1 or more")
  if(pedtools::is.ped(pm)){
    labelspm = labels(pm)
    if(!all(vp %in% labels(pm)))
      errortext = c("Third argument must be a subset of labels of first argument")
  }
  if(pedtools::is.ped(am)){
    labelsam = labels(am)
    if(!all(mp %in% labels(am)))
      errortext = c("Fourth argument must be a subset of labels of second argument")
  }
  if(pedtools::is.pedList(pm)){
    labelspm = unlist(lapply(pm, labels))
    if(!all(vp %in% labelspm))
      errortext = c("Third argument must be a subset of labels of first argument")
  }
  if(pedtools::is.pedList(am)){
    labelsam = unlist(lapply(am, labels))
    if(!all(mp %in% labelsam))
      errortext = c("Fourth argument must be a subset of labels of second argument")
  }
  if(length(unique(vp)) < nvp )
    errortext = c("Victims must be unique")
  if(length(unique(mp)) < nmp )
    errortext = c("Missing persons must be unique")
  g = getAlleles(pm)
  if(all(is.na(g)))
    errortext = c("No genotypes for victims")
  g = getAlleles(am, ids = mp)
  if(!all(is.na(g)))
    errortext = c("Genotypes for missing persons")
  likNULL = prod(LR(list(pm, am), 1)$likelihoodsPerSystem)
  if (likNULL == 0){
    errortext = "Initial data has 0 likelihood"
   }
  list(error = errortext, lik0 = likNULL)
}