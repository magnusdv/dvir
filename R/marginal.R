#' Check and evaluate individual candidates
#' 
#' Each potential move is evaluated and the search space reduced. 
#' 
#' @param from A list of singletons, the victims. 
#' @param to A list of pedigrees. The reference families.
#' @param ids.to Character vector with names of missing persons.
#' @param moves List with possible marginal moves.
#' @param limit Double. Lower threshold for LR.
#' @param verbose Logical.
#' @param sorter Logical, sorts output according to LR.
#' @param nkeep integer. No of moves to keep, all if `NULL`.
#' @details The potential reduction only affects the list of moves returned, all LRs are kept.
#' Specifying `nkeep`can give further reduction.
#' 
#' 
#' @return A list with moves and log likelihoods.
#' @import forrel
#' @importFrom pedprobr likelihood
#' @export
#' @examples
#' \dontrun{
#' library(forrel)
#' n = 7
#' ids.from = paste("V", 1:n, sep = "")
#' sex = c(rep(1, n-1), 2)
#' df = data.frame(id = ids.from, fid = 0, mid = 0, sex = sex, 
#'                 a1 = c(1,2,1,1,2,2,2), 
#'                 a2 = c(1,2,1,1,2,2,2))
#' locus_annotations = list(alleles = 1:3, afreq = c(1, 1, 1)/3)
#' from = as.ped(df, locusAttributes = locus_annotations)
#' names(from) = ids.from
#' to = nuclearPed(3, father = "R1", mother = "R2", children= c("MP1","MP2","MP3"))
#' m = marker(to, "R1" = 1, "R2" = 1,   alleles = 1:3, afreq = c(1, 1, 1)/3, name = "a1")
#' to = addMarkers(to, m)
#' ids.to = c("MP1", "MP2","MP3")
#' plotPedList(list(from, to), marker = 1)
#' limit = -1; nbest = NULL; extend = T; merge = FALSE; rename = TRUE;moves = NULL

#' moves = list(V1 = c("MP1", "V1", "MP2"), V3 = c("MP1","MP2", "MP3"), V4 = c("V4", "MP3"), 
#'              V7 = c("V7"))
#' # all 4*6+1 possible moves ignoring sex
#' moves = list(V1 = c("V1", ids.to), V2 = c("V2", ids.to), V3 = c("V3", ids.to),
#'              V4 = c("V4", ids.to), V5 = c("V5", ids.to), V6 = c("V6", ids.to),
#'              V7 = c("V7"))
#'             
#' res = marginal(from, to,  ids.to, moves, limit = 0, verbose = T, sorter = T, nkeep= 2)
#' res = marginal(from, to,  ids.to, moves, limit = -1, sorter = T,  nkeep = 3)
#' res2 = global(from, to, ids.to, moves = res[[1]], limit = 0)
#' moves = list(V1 = c("V1", "MP1", "MP2"))
#' res = marginal(from, to,  ids.to, moves, limit = 1)
#' }

marginal = function(from, to, ids.to, moves, limit = 0.1, 
                    verbose = FALSE, sorter = FALSE, nkeep = NULL){
  if( !sorter & !is.null(nkeep))
    stop("If keep is not NULL, sorter should be TRUE")
  if(is.null(moves)) # Generate moves
    moves = generateMoves(from = from, to = to,  ids.to = ids.to)
  else # remove elements with missing
    moves = moves[unlist(lapply(moves, function(x) !any(is.na(x))))]
  res = checkDVI(from = from, to = to,  ids.to = ids.to , moves = moves)

  ids.from = as.character(lapply(from, function(x) x$ID))
  marks = 1:nMarkers(from)
  names(from) = ids.from
  loglik0 = sum(likelihood(from, marker = marks, logbase = exp(1)), eliminate = 1) +
            sum(likelihood(to,   marker = marks, logbase = exp(1)), eliminate = 1)
  LR = list()
  if(loglik0 == -Inf)
    stop("Impossible initial data")
  LR = list()
  n.moves = length(moves)
  res = list()
  moves.reduced = moves
  for (i in 1:n.moves){
    if(verbose) cat("Move: ",i, "of", n.moves, "\n")
    res[[i]] = screen1(from, to, ids.to, moves, vict = i, 
                       loglik0 = loglik0, LRlimit = limit, verbose = verbose,
                       sorter = sorter, nkeep = nkeep)
    moves.reduced[[i]] = res[[i]]$move1Kept
    LR[[i]] = res[[i]]$LR
  }
  list(moves = moves.reduced, LR = LR)
}
    
                  
screen1 = function(from, to, ids.to, moves, loglik0, vict = 1, LRlimit = 0.1, 
                   verbose = F, sorter = FALSE, nkeep = NULL){
  marks = 1:nMarkers(from)
  idFrom = names(moves[vict])  # id of victim to potentially move
  from1 = from[[idFrom]] # singleton corresponding to idFrom
  move1 = moves[[vict]] # ids of potential moves
  idTo = as.character(move1)
  nm = length(move1)
  logl = keep = LR = rep(NA, nm)
  for (i in 1:nm){
    if(verbose) cat("Iteration ",i, "of", nm, "\n")
    idTo = move1[i]
    if(idTo == idFrom){ 
        logl[i] = loglik0
        keep[i] = LRlimit <1
        LR[i] = 1
    } else {
        from2 = relabel(from1, idTo, idFrom)
        am2 = transferMarkers(from2, to, idTo, erase = FALSE)
        pm2 =  setAlleles(from, ids = idFrom, alleles = 0)
        logl[i] = sum(likelihood(pm2, marker = marks, logbase = exp(1)), eliminate = 1) +
                  sum(likelihood(am2, marker = marks, logbase = exp(1)), eliminate = 1)
        LR[i] = exp(logl[i]-loglik0) 
        keep[i] = LR[i] > LRlimit
      }
  }
  names(logl) = names(LR) = names(keep) = move1
  if(sorter){
   ord = order(LR, decreasing = TRUE)
   logl = logl[ord]
   LR = LR[ord]
   keep = keep[ord]
   move1 = move1[ord]
   if(!is.null(nkeep)){
     nk = min(nkeep, length(ord))
     keepMoves = (1:nk)[keep[1:nk]]
     move1 = move1[keepMoves]
   }  else
    move1 = move1[keep]
  }
    if (length(move1) == 0){
      if(verbose) cat("move not possible for: ", idFrom,"\n")
      move1 = NA
    }
    list(loglik = logl, LR = LR, keep = keep, 
              move1Kept = move1)
}
  


