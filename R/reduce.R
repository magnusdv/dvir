#' Finds persons that can be removed
#' 
#' LR is calculated for each possible move without updating.
#' Each line of the output table gives the result for 
#' moves meeting the LR requirement. The lists of relevant 
#' victims and missing are returned. 
#' 
#' 
#' @param from A list of singeltons 
#' @param to A list of pedigrees
#' @param ids.from Character vector with names of victims
#' @param ids.to Character vector with names of missing persons
#' @param likNull Double
#' @param limit Double
#' @return A list of three components
#' 
#' * eliminateTable Data frame. A line for each move with LR above limit
#' 
#' * ids.from Character vector. id-s of samples in ids.from not eliminated
#' 
#' * ids.to Character vector. id-s of samples in ids.to not eliminated

#' @export
#' @examples 
#' 
#' als = 1:2
#' freq = c(0.5, 0.5)
#' mutmodel = "proportional"
#' rate = 0.000
#' v1 = singleton("V1")
#' m = marker(v1, alleles = als, afreq = freq, "V1" = c(1, 1),
#' mutmod = mutmodel, rate = rate)
#' v1 = addMarkers(v1, m)
#' v2 = singleton("V2")
#' m = marker(v2, alleles = als, afreq = freq, "V2" = c(1, 2),
#' mutmod = mutmodel, rate = rate)
#' v2 = addMarkers(v2, m)
#' v3 = singleton("V3", 2)
#' m = marker(v3, alleles = als, afreq = freq, "V3" = c(2, 2),
#'           mutmod = mutmodel, rate = rate)
#' v3 = addMarkers(v3, m)
#' from = list(v1, v2, v3)
#' to = nuclearPed(3, sex = c(1,1,2), father = "fa", mother = "mo", children = c("MP1", "MP2", "MP3"))
#' m = marker(to, alleles = als, afreq = freq, "fa" = 1:2, "mo" = c(1,1),
#'         mutmod = mutmodel, rate = rate)
#' to = addMarkers(to, m)
#' plotPedList(list(from, to), marker = 1)
#' ids.from = c("V1", "V2", "V3")
#' ids.to = c("MP1","MP2", "MP3")
#' check = checkInput(from, to, ids.from, ids.to)
#' if(!is.null(check$error)) stop(check$error)
#' lik0 = check$lik0
#' limit = 0
#' res = reduce(from, to, ids.from, ids.to,likNULL = lik0, limit = limit)
#' 

#' pm1 = singleton("V1")
#' m = marker(pm1, "V1" = 1, alleles = 1:2)
#' pm1 = addMarkers(pm1, m)
#' pm2 = singleton("V2")
#' m = marker(pm2, "V2" = 2, alleles = 1:2)
#' pm2 = addMarkers(pm2, m)
#' from = list(pm1, pm2)
#' ids.from = c("V1", "V2")
#' to = nuclearPed(1, father = "MP1", mother = "MP2", children= "R1")
#' m = marker(to, "R1" = 1, alleles = 1:2)
#' to = addMarkers(to, m)
#' ids.to =c("MP1", "MP2")
#' plotPedList(list(from, to), marker = 1)
#' check = checkInput(from, to, ids.from, ids.to)
#' lik0 = check$lik0
#' limit = 0
#' res = reduce(from, to, ids.from, ids.to,likNULL = lik0, limit = limit)

reduce = function(from, to, ids.from, ids.to, likNULL = 1, limit = 0){
  tab = generate(from, to, ids.from, ids.to)$onestep
  if(is.null(tab))
    return(NULL)
  liks = apply(tab ,1, function(x, from, to) 
    lik1(from, to, x[1] , x[2])$lik, from, to)
  lr = liks/likNULL
  index = lr > limit
  ni = sum(index)
  if (ni > 0 ){
    tab2 = matrix(tab[index,], nrow = ni)
    ret = matrix(cbind(tab2, liks[index], rep(likNULL, ni), lr[index]), ncol = 5)
    colnames(ret) = c("from", "to", "lik", "lik0", "LR")
    ids.from = sort(unique(ret[,1]))
    ids.to = sort(unique(ret[,2]))
    ret = list(eliminateTable = data.frame(ret), 
               ids.from = ids.from, ids.to = ids.to)
  } else
    return(NULL)
  if (length(ids.from) < 1 | length(ids.to) < 1)
    return(NULL)
  ret
}