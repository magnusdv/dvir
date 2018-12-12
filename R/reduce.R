#' remove irelevant victims and missing persons
#' 
#' @param pm A list of singeltons 
#' @param am A list of pedigrees
#' @param vp Character vector with names of victims
#' @param mp Character vector with names of missing persons
#' @param likNull Double
#' @param limit Double
#' @export
#' @examples 
#' pm1 = singleton("V1")
#' m = marker(pm1, "V1" = 1, alleles = 1:2)
#' pm1 = addMarkers(pm1, m)
#' pm2 = singleton("V2")
#' m = marker(pm2, "V2" = 2, alleles = 1:2)
#' pm2 = addMarkers(pm2, m)
#' p = list(pm1, pm2)
#' vp = c("V1", "V2", "V3", "V4")
#' am = nuclearPed(1, father = "MP1", mother = "MP2", children= "R1")
#' m = marker(am, "R1" = 1, alleles = 1:2)
#' am = addMarkers(am, m)
#' mp =c("MP1", "MP2")
#' reduce(pm, am, vp, mp)
#' 

reduce = function(pm, am, vp, mp, likNULL = 1, limit = 0){
  tab = generate(pm, am, vp, mp)$onestep
  if(is.null(tab))
    return(NULL)
  liks = apply(tab ,1,function(x, pm, am) 
    lik1(x[1] , x[2], pm, am)$lik, pm, am)
  lr = liks/likNULL
  index = lr > limit
  ni = sum(index)
  if (ni > 0 ){
    ret = matrix(cbind(tab[index, ], liks[index], rep(likNULL, ni), lr[index]), ncol = 5)
    colnames(ret) = c("from", "to", "lik", "lik0", "LR")
    vp = sort(unique(ret[,1]))
    mp =sort(unique(ret[,2]))
    ret = list(ret = data.frame(ret), vp = vp, mp =mp)
  } else
    return(NULL)
  if (length(vp) < 1 | length(mp) < 1)
    return(NULL)
  ret
}