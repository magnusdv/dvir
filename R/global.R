#' Search for best DVI solution
#' 
#' Victims are given as a list of singletons and references as list of pedigrees.
#' Based on a specification of moves for each victim, if any, all possible moves 
#' are evaluated and a solution ranked according to the likelihood is returned. 
#' 
#' @param from A list of singletons. 
#' @param to A list of pedigrees.
#' @param ids.to Character vector with names of missing persons.
#' @param moves List of length equal from with possible marginal moves.
#' @param limit Double. Only moves with LR above limit are kept.
#' @param verbose Logical.
#' @details This is currently a brute force appraoch, all possibilities are evaluated
#' 
#' @return A data frame. Each row describes an a priori possible move. The log likelihood and 
#' the null likelihood is given as well as the LR with the null hypothesis in the numerator.
#' There is limited checking. For checking, preprocess with [marginal()]
#' @seealso [marginal()]
#' @export
#' @examples
#' \dontrun{
#' library(forrel)
#' library(pedprobr)
#' n = 7
#' ids.from = paste("V", 1:n, sep = "")
#' sex = c(rep(1, n-1), 2)
#' df = data.frame(id = ids.from, fid = 0, mid = 0, sex = sex, 
#'                 a1 = c(1,2,1,1,2,2,2), 
#'                 a2 = c(1,2,1,1,2,2,2))
#' locus_annotations = list(alleles = 1:3, afreq = c(1, 1, 1)/3)
#' from = as.ped(df, locusAttributes = locus_annotations)
#' to = nuclearPed(3, father = "R1", mother = "R2", children= c("MP1","MP2","MP3"))
#' m = marker(to, "R1" = 1, "R2" = 1,   alleles = 1:3, afreq = c(1, 1, 1)/3, name = "a1")
#' to = addMarkers(to, m)
#' ids.to = c("MP1", "MP2","MP3")
#' plotPedList(list(from, to), marker = 1)
#' moves = list(V1 = c("V1", "MP1", "MP2"), V3 = c("MP1","MP2", "MP3"), V4 = c("V3", "MP3"), 
#'              V7 = c("V7","MP3"))
#' res = global(from, to, ids.to, moves, limit = 0)
#' moves = moves[1:2]
#' res = global(from, to, ids.to, moves, limit = 0)
#' moves = moves[1]
#' res = global(from, to, ids.to, moves, limit = 0)
#' 
#' # Example 2
#' # Checked against; http://familias.name/BookKEP/globalExample2.fam
#' Several reference families and relationship between victims
#' V1 = singleton("V1")
#' m = marker(name = "m1", V1, alleles = 1:10, "V1" = c(1,1))
#' V1 = setMarkers(V1, m)
#' V2 = singleton("V2")
#' m = marker(name = "m1", V2, alleles = 1:10, "V2" = c(1,2))
#' V2 = setMarkers(V2, m) 
#' V3 = singleton("V3",2)
#' m = marker(name = "m1", V3, alleles = 1:10, "V3" = c(3,4))
#' V3 = setMarkers(V3, m) 
#' V4 = singleton("V4")
#' m = marker(name = "m1", V4, alleles = 1:10, "V4" = c(3,4))
#' V4 = setMarkers(V4, m)  
#' from = list(V1 = V1, V2 = V2, V3 = V3, V4 = V4)
#' AM1 = nuclearPed(1, father = "MP1", child = "MP2", mother ="R1")
#' m = marker(name = "m1",AM1, alleles = 1:10, "R1" = c(2, 2))
#' AM1 = setMarkers(AM1, m)
#' AM2 = halfSibPed(sex1 = 1, sex2 = 2)
#' AM2 = relabel(AM2, c("R2", "MO2", "FA2", "MP3", "MP4"), 1:5)
#' m = marker(name = "m1", AM2, alleles = 1:10, "R2" = c(3,3))
#' AM2 = setMarkers(AM2, m)
#' to = list(AM1 = AM1, AM2 = AM2)
#' plotPedList(list(from, to), marker = "m1", shaded = typedMembers)
#' ids.to = c("MP1", "MP2", "MP3", "MP4") # user specified
#' males = ids.to[1:3]
#' females = ids.to[4]
#' moves = list(V1 = c("V1", males), 
#'              V2 = c("V2", males),
#'              V3 = c("V3", females), 
#'              V4 = c("V4", males))
#' res = global(from, to, ids.to, moves = moves, limit = 0, verbose = T)
#' resmarg = marginal(from, to, ids.to, moves = moves, limit = 0, 
#'              verbose = T, sorter = TRUE, nkeep = 2)
#' res = global(from, to, ids.to, moves = resmarg[[1]], limit = 0, verbose = T)
#' 
#' #From fam file
#' x = readFam("http://familias.name/BookKEP/globalExample2.fam")
#' from = x[[1]][1:4]
#' to = x[[1]][5:6]
#' plotPedList(list(from,to))
#' ids.to = c("MP1", "MP2", "MP3", "MP4") # user specified
#' males = ids.to[1:3]
#' females = ids.to[4]
#' moves = list(V1 = c("V1", males), V2 = c("V2", males),
#'              V3 = c("V3", females), V4 = c("V4", males))
#' res = global(from, to, ids.to, moves = moves, limit = 0, verbose = T)
#' }

global = function(from, to,  ids.to, moves, limit = 0, verbose = F){
  if(length(ids.to) < 1)
    stop("Third vector must a vector with names of missing persons")
  if(!is.list(moves))
    stop("Fourth argument must be a list")
  ids.from = as.character(lapply(from, function(x) x$ID))
  # names(from) = ids.from # Names of victims, singletons
  marks = 1:nMarkers(from)
  loglik0 = sum(likelihood(from, marker = marks, logbase = exp(1))) +
            sum(likelihood(to,   marker = marks, logbase = exp(1)))
  if(loglik0 == -Inf)
    stop("Impossible initial data")
  tomove = names(moves) #names of all victims up for a potential move
  moves2 = expand.grid.nodup(moves) #each element of moves2 is a possible move
  nm = length(moves2)
  if(nm == 0)
    stop("No possible solutions specified, possibly identical moves")
  
  loglik = LR = rep(NA, nm)

  for (i in 1:nm){
    if(verbose) cat("iteration ", i, "of", nm," ","\n")
    thisMove = setdiff(moves2[[i]], tomove)
    if(length(thisMove) == 0){
      loglik[i] = loglik0
      moves2[[i]] = data.frame(moves2[[i]], loglik = loglik[i],loglik0 = loglik0, LR = 1)
    } else {
      idFrom = names(thisMove)
      from2 = relabel(from, thisMove[1,], names(thisMove))
      to2 = transferMarkers(from2, to, erase  = FALSE)
      from2 = setAlleles(from, ids = idFrom, alleles = 0)
      loglik[i] = sum(likelihood(from2, marker = marks, logbase = exp(1))) +
                  sum(likelihood(to2,   marker = marks, logbase = exp(1)))
      moves2[[i]] = data.frame(moves2[[i]], loglik = loglik[[i]], 
                              loglik0 = loglik0, LR = exp(loglik[i] - loglik0))
    }
  }
  ord = order(loglik, decreasing = T)
  moves2 = moves2[ord]
  keep = unlist(lapply(moves2, function(x) x$LR > limit))
  moves2 = moves2[keep]
  moves2 = data.frame(t(sapply(moves2,c)))
  colnames(moves2) = c(tomove, c("loglik", "loglik0", "LR"))
  moves2
}

