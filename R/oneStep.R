#' Forward DVI selection, one at a time
#' 
#' @param from A list of singeltons 
#' @param to A list of pedigrees
#' @param ids.from Character vector with names of victims
#' @param ids.to Character vector with names of missing persons
#' @param LRlimit Double. Threshold for LR
#' @param lik0 Double. 
#' @return output
#' @examples
#' ids.from = c("V1", "V2")
#' from = singletonList(2, ids = ids.from)
#' m = marker(from[[1]], alleles = 1:2, "V1" = 2)
#' from[[1]] = addMarkers(from[[1]], m)
#' m = marker(from[[2]], alleles = 1:2, "V2" = 1:2)
#' from[[2]] = addMarkers(from[[2]], m)
#' ids.to = c("MP1", "MP2")
#' to = nuclearPed(2, children = ids.to, father = "R1", mother = "R2")
#' m = marker(to, "R1" = 1, "R2"= 1:2, alleles = 1:2)
#' to = addMarkers(to, m)
#' plotPedList(list(from, to), marker = 1)
#' lik0 = prod(LR(list(from, to), 1)$likelihoodsPerSystem)
#' limit = -1
#' one = oneStep(from, to, ids.from, ids.to, LRlimit = limit, lik0 = lik0)
#' @export
oneStep <-
  function(from, to, ids.from , ids.to, LRlimit = 1, lik0 = NULL, 
           maxid = NULL){
  maxno = min(length(ids.from), length(ids.to))  
  mi = ifelse(is.null(maxid), maxno, min(maxid, maxno))
  # Initalise
  keep1 = matrix(nrow = mi+1, ncol = 2)
  colnames(keep1) = c("from", "to")
  keep2 = matrix(nrow = mi+1, ncol = 4)
  colnames(keep2) = c("lik", "LR", "prior", "posterior")
  keep2[1,1] = lik0
  i = 1
  lr = LRlimit + 1
  tab = generate(from, to, ids.from, ids.to)$onestep
  if(is.null(tab)) return(NULL)
  tab = matrix(generate(from, to, ids.from, ids.to)$onestep, ncol = 2)
  while(i <= mi & lr > LRlimit & !is.null(tab)){
      res1 = apply(tab, 1, function(x, from, to) 
        lik1(from, to, x[1] , x[2]), from, to)
      likres = unlist(lapply(res1, function(x) x$lik))
      no = which.max(likres)
      likMax = likres[no]
      res1 = res1[[no]]
      from1 = as.character(tab[no, 1])
      to1 = as.character(tab[no, 2])
      
      ids.from = setdiff(ids.from, from1)
      ids.to = setdiff(ids.to,to1)
      keep1[i+1,1] = from1
      keep1[i+1,2] = to1
      keep2[i+1,1] = likMax
      lr = likMax/keep2[i,1] 
      keep2[i+1,2] = lr
      keep2[i+1,3] = NA
      keep2[i+1,4] = likMax/sum(likres)
      index = !tab[,1] %in% from1 & !tab[,2] %in% to1
      if(sum(index) == 0)
        tab = NULL
      else
        tab = matrix(tab[index, ], ncol = 2)
      if(is.na(lr)) lr = LRlimit - 1
      if(LRlimit < 0 | lr > LRlimit) {
        from = res1$pm
        to = res1$am
        i = i + 1
      }
  }
  res = data.frame(keep1, keep2)[1:i, ]
  rownames(res) = paste("step", 0:(i-1), sep ="")
  list(singleStep = res, pm2 = from, am2 = to)
  }
