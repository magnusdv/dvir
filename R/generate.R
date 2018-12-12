#' Find next candidates
#' 
#' New moves are generated, NULL
#' if none are possible, for one and two step
#' 
#' @param pm A list of singeltons 
#' @param am A list of pedigrees 
#' @param vp Character vector with names of victims
#' @param mp Character vector with names of missing persons
#' @param one logic If TRUE, only one-step
#' @return list of matrices of one and two-moves
#' @examples 
#' vp = c("V1", "V2")
#' pm = singletonList(2, ids = vp)
#' m = marker(pm[[1]], alleles = 1:2, "V1" = 1)
#' pm[[1]] = addMarkers(pm[[1]], m)
#' m = marker(pm[[2]], alleles = 1:2, "V2" = 1:2)
#' pm[[2]] = addMarkers(pm[[2]], m)
#' mp = c("MP1", "MP2")
#' am = nuclearPed(2, children = mp, father = "R1", mother = "R2")
#' m = marker(am, "R1" = 1:2, "R2"= 1:2)
#' am = addMarkers(am, m)
#' # The thre possibilities: a move of two, one or none:
#' tab = generate(pm, am, vp, mp)
#' pm[[1]] = swapSex(pm[[1]], "V1")
#' tab = generate(pm, am, vp,mp)
#' pm[[2]] = swapSex(pm[[2]], "V2")
#' tab = generate(pm, am, vp,mp)
#' @export
#' 
generate = function(pm, am, vp, mp, one = TRUE){
  nvp = length(vp)
  nmp = length(mp)
  if(one)
    tab2  = NULL
  else
    tab2 = loopPair(vp, mp, pm, am)
  #tab2 = ifelse(one, NULL, loopPair(vp, mp, pm, am))
  amsex = data.frame(mp = mp, sex = getSex(am, mp))
  rownames(amsex) = mp
  pmsex = data.frame(vp = vp, sex = getSex(pm, vp))
  rownames(pmsex) = vp
  sex = list(pmsex = pmsex, amsex = amsex)
  tab1 = loop2(vp,mp, sex)
  list(onestep = tab1, twostep = tab2)
}
