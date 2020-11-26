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
#' @param nkeep integer. No of moves to keep, all if `NULL`.
#' @details The potential reduction only affects the list of moves returned, all LRs are kept.
#' Specifying `nkeep`can give further reduction.
#' 
#' 
#' @return A list with moves and log likelihoods.
#' 
#' @importFrom pedprobr likelihood
#' @export
#' @examples
#' \donttest{
#' library(pedtools)
#' 
#' # Attributes of a single marker
#' loc = list(name = "m", alleles = 1:3)
#' 
#' # Victims
#' vics = paste0("V", 1:7)
#' sex = c(1, 1, 1, 1, 1, 1, 2)
#' df = data.frame(id = vics, fid = 0, mid = 0, sex = sex, 
#'                 m = c("1/1", "2/2", "1/1", "1/1", "2/2", "2/2", "2/2"))
#' from = as.ped(df, locusAttributes = loc)
#' names(from) = vics
#' 
#' # Reference families
#' MPs = c("MP1", "MP2", "MP3")
#' to = nuclearPed(3, father = "R1", mother = "R2", children = MPs)
#' data = data.frame(m = c("1/1", "1/1"), row.names = c("R1", "R2"))
#' to = setMarkers(to, alleleMatrix = data, locusAttributes = loc)
#' 
#' plotPedList(list(from, to), marker = 1)
#' 
#' # moves = list(V1 = c("MP1", "V1", "MP2"), 
#' #              V3 = c("MP1", "MP2", "MP3"), 
#' #              V4 = c("V4", "MP3"), 
#' #              V7 = c("V7"))
#' #             
#' # # all 4 * 6 + 1 =possible moves ignoring sex
#' # moves = list(V1 = c("V1", MPs), V2 = c("V2", MPs), V3 = c("V3", MPs),
#' #              V4 = c("V4", MPs), V5 = c("V5", MPs), V6 = c("V6", MPs),
#' #              V7 = c("V7"))
#' #
#' # res = marginal(from, to, MPs, moves, limit = 0, verbose = TRUE, nkeep= 2)
#'  
#' moves = generateMoves(from, to, MPs)
#' 
#' res = marginal(from, to, MPs, moves, limit = -1, nkeep = 3)
#' res2 = global(from, to, MPs, moves = res[[1]], limit = 0)
#' # moves = list(V1 = c("V1", "MP1", "MP2"))
#' # res = marginal(from, to,  MPs, moves, limit = 1)
#' }

marginal = function(from, to, ids.to, moves, limit = 0.1, 
                    verbose = FALSE, nkeep = NULL){
  
  if(is.null(moves)) # Generate moves
    moves = generateMoves(from = from, to = to,  ids.to = ids.to)
  
  # Check consistency
  res = checkDVI(from = from, to = to, ids.to = ids.to, moves = moves)

  marks = 1:nMarkers(from)
  
  # Ensure correct names
  names(from) = unlist(labels(from), use.names = FALSE)
  vics = names(moves) # normally in the same order as names(from)
  
  # Loglik of each victim
  logliks.PM = vapply(from, loglikTotal, markers = marks, FUN.VALUE = 1)
  
  # log-likelihood of H0
  loglik0 = sum(logliks.PM) + loglikTotal(to, marks)
  
  if(loglik0 == -Inf)
    stop("Impossible initial data")
  
  # For each victim, compute marginal LRs of each move
  LR.list = lapply(vics, function(v) {
    
    # Vector of moves for v
    mps = moves[[v]]
    
    # Corresponding vector of LRs
    lrs = vapply(moves[[v]], function(mp) {
      if(mp == "*") 
        return(1)
      
      # Likelihood of remaining PMs
      loglik.remaining = sum(logliks.PM[setdiff(vics, v)])
      
      # Likelihood of families after move
      am2 = transferMarkers(from[[v]], to, idsFrom = v, idsTo = mp, erase = FALSE)
      loglik.fam = loglikTotal(am2, marks)
      
      # Total loglik after move
      loglik.move = loglik.remaining + loglik.fam
      
      # Return LR of move
      exp(loglik.move - loglik0)
    }, FUN.VALUE = numeric(1))
    
    # Return sorted vector
    sort(lrs, decreasing = TRUE)
  })
  
  names(LR.list) = vics
  
  # Matrix of marginal LRs (filled with 0's)
  LR.table = matrix(0, nrow = length(vics), ncol = length(ids.to), 
                    dimnames = list(vics, ids.to))
  
  # Fill matrix row-wise
  for (v in vics) {
    lrs = LR.list[[v]]
    lrs = lrs[names(lrs) != "*"]  # remove do-nothing move
    LR.table[v, names(lrs)] = unname(lrs)
  }
  
  # Reduce moves according to `limit` and/or nkeep
  moves.reduced = lapply(LR.list, function(lrs) {
    newmoves = names(lrs)[lrs >= limit]
    if(!is.null(nkeep) && length(newmoves) > nkeep)
      length(newmoves) = nkeep
    newmoves
  })
  
  list(moves = moves.reduced, LR.list = LR.list, LR.table = LR.table)
}
   


# Not used 
screen1 = function(from, to, ids.to, moves, loglik0, vict = 1, LRlimit = 0.1, 
                   verbose = F, sorter = FALSE, nkeep = NULL){
  ids.from = unlist(labels(from))
  marks = 1:nMarkers(from)
  idFrom = names(moves[vict])  # id of victim to potentially move
  from1 = from[[idFrom]] # singleton corresponding to idFrom
  move1 = moves[[vict]] # ids of potential moves
  idTo = as.character(move1)
  nm = length(move1)
  
  # Loglik of each victim
  logliks.PM = vapply(from, loglikTotal, markers = marks, FUN.VALUE = 1)
  names(logliks.PM) = ids.from
  
  logl = keep = LR = rep(NA, nm)
  for (i in 1:nm){
    if(verbose) cat("Iteration ",i, "of", nm, "\n")
    idTo = move1[i]
    if(idTo == idFrom){ 
        logl[i] = loglik0
        keep[i] = LRlimit <1
        LR[i] = 1
    } else {
      # Likelihood of remaining PMs
      ids.remaining = setdiff(ids.from, idFrom)
      loglik.remaining = sum(logliks.PM[ids.remaining])
      
      # Likelihood of families after move
      am2 = transferMarkers(from, to, idsFrom = idFrom, 
                            idsTo = idTo, erase = FALSE)
      loglik.fam = loglikTotal(am2, marks)
      
      logl[i] = loglik.remaining + loglik.fam
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
  


