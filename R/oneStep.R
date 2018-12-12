#' Forward DVI selection, one at a time
#' 
#' @param pm A list of singeltons 
#' @param am A list of pedigrees
#' @param vp Character vector with names of victims
#' @param mp Character vector with names of missing persons
#' @param LRlimit Double. Threshold for LR
#' @param lik0 Double. 
#' @return output
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
#' lik0 = prod(LR(list(pm, am), 1)$likelihoodsPerSystem)
#' limit = -1
#' one = oneStep(pm, am, vp, mp, LRlimit = limit, lik0 = lik0)
#' @export
oneStep <-
  function(pm, am, vp , mp, LRlimit = 1, lik0 = NULL, maxid = NULL){
  maxno = min(length(vp), length(mp))  
  if(is.null(maxid))
    mi = maxno
  else
    mi = min(maxid, maxno)
  # Initalise
  keep1 = matrix(nrow = mi+1, ncol = 2)
  colnames(keep1) = c("from", "to")
  keep2 = matrix(nrow = mi+1, ncol = 4)
  colnames(keep2) = c("lik", "LR", "prior", "posterior")
  keep2[1,1] = lik0
  i = 1
  lr = LRlimit + 1
  tab = generate(pm, am, vp, mp)$onestep
  if(!is.null(tab))
    tab = matrix(generate(pm, am, vp, mp)$onestep, ncol = 2)
  while(i <= mi & lr > LRlimit &! is.null(tab)){
      res1 = apply(tab, 1, function(x, pm, am) 
        lik1(x[1] , x[2], pm, am), pm, am)
      likres = unlist(lapply(res1, function(x) x$lik))
      no = which.max(likres)
      likMax = likres[no]
      res1 = res1[[no]]
      from = as.character(tab[no, 1])
      to = as.character(tab[no, 2])
      pm = res1$pm
      am = res1$am
      vp = setdiff(vp, from)
      mp = setdiff(mp,to)
      keep1[i+1,1] = from
      keep1[i+1,2] = to
      keep2[i+1,1] = likMax
      lr = likMax/keep2[i,1]
      keep2[i+1,2] = lr
      keep2[i+1,3] = NA
      keep2[i+1,4] = likMax/sum(likres)
      index = !tab[,1] %in% from & !tab[,2] %in% to
      if(sum(index) == 0)
        tab = NULL
      else
        tab = matrix(tab[index, ], ncol = 2)
      if(is.na(lr)) lr = LRlimit - 1
      if(lr > LRlimit) i = i + 1
  }
  res = data.frame(keep1, keep2)[1:i, ]
  rownames(res) = paste("step", 0:(i-1), sep ="")
  list(result = res, pm2 = pm, am2 = am)
  }