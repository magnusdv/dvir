#' Move individuals, find likelihoods and update
#' 
#' Individuals are moved from pm to am data.
#' The likelihood is calculated and pm and am
#' are updated 
#' 
#' @param from Character vector with names of victims
#' @param to Character vector with names of missing persons
#' @param pm A list of singeltons 
#' @param am A list of pedigrees 
#' @return Likelihood and updated pm and am lists
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
#' res = lik1(vp, mp, pm, am) 
#' res = lik1(vp[1], mp[1], pm, am) 
#' @export
lik1 <-
function(from, to, pm, am){
    g = getAlleles(pm, ids = from)
    pm2 = setAlleles(pm, ids = from, alleles = 0 )
    rownames(g) = to
    am2 = setAlleles(am, id  = to, alleles = g)
    list(lik = prod(LR(list(pm2, am2),1)$likelihoodsPerSystem),
         pm2 = pm2, am2 = am2)
    }
