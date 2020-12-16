#' Single-search LR matrix
#'
#' For a given DVI problem, compute the LR matrix consisting of individual
#' likelihood ratios \eqn{LR_{i,j}} comparing the assignment \eqn{V_i = M_j} to
#' the null. The output may be reduced by specifying arguments `limit` or
#' `nkeep`.
#'
#' @param pm A list of singletons, the victims.
#' @param am A list of pedigrees. The reference families.
#' @param missing A character vector with names of missing persons.
#' @param moves A list with possible assignments.
#' @param limit A positive number: only single-search LR values above this are
#'   considered.
#' @param nkeep An integer. No of moves to keep, all if `NULL`.
#' @param check A logical, indicating if the input data should be checked for
#'   consistency.
#' @param verbose A logical.
#'
#' @return A list with moves and log likelihoods.
#'
#' @examples
#'
#' pm = example1$pm
#' am = example1$am
#' missing = example1$missing
#' 
#' singleSearch(pm, am, missing)
#'
#' @export
singleSearch = function(pm, am, missing, moves = NULL, limit = 0, nkeep = NULL, 
                    check = TRUE, verbose = FALSE){
  
  if(is.singleton(pm)) pm = list(pm)
  if(is.ped(am)) am = list(am)
  
  if(is.null(moves)) # Generate moves
    moves = generateMoves(pm, am, missing = missing, expand.grid = FALSE)
  
  # Check consistency
  if(check)
    checkDVI(pm, am, missing = missing, moves = moves)

  # Ensure correct names
  names(pm) = unlist(labels(pm), use.names = FALSE)
  vics = names(moves) # normally in the same order as names(pm)
  
  marks = 1:nMarkers(pm)
  
  # Loglik of each victim
  logliks.PM = vapply(pm, loglikTotal, markers = marks, FUN.VALUE = 1)
  
  # Loglik of each ref family
  logliks.AM = vapply(am, loglikTotal, markers = marks, FUN.VALUE = 1)
  
  # log-likelihood of H0
  loglik0 = sum(logliks.PM) + sum(logliks.AM)
  
  if(loglik0 == -Inf)
    stop("Impossible initial data")
  
  # For each victim, compute the LR of each move
  LR.list = lapply(vics, function(v) {
    
    # Vector of moves for v
    missing = moves[[v]]
    
    # Corresponding vector of LRs
    lrs = vapply(moves[[v]], function(mp) {
      if(mp == "*") 
        return(1)
      
      # Make copy of AM likelihoods (vector)
      logliks.AM.new = logliks.AM
      
      # The relevant AM component 
      compNo = getComponent(am, mp, checkUnique = TRUE)
      
      # Move victim data to `mp`
      comp = transferMarkers(pm[[v]], am[[compNo]], idsFrom = v, idsTo = mp, erase = FALSE)
      
      # Update likelihood of this comp
      logliks.AM.new[compNo] = loglikTotal(comp, marks)
      
      # Likelihood of remaining PMs
      logliks.PM.new = logliks.PM[setdiff(vics, v)]
      
      # Total loglik after move
      loglik.move = sum(logliks.PM.new) + sum(logliks.AM.new)
      
      # Return LR of move
      exp(loglik.move - loglik0)
    }, FUN.VALUE = numeric(1))
    
    # Return sorted vector
    sort(lrs, decreasing = TRUE)
  })
  
  names(LR.list) = vics
  
  # Matrix of individual LRs (filled with 0's)
  LR.table = matrix(0, nrow = length(vics), ncol = length(missing), 
                    dimnames = list(vics, missing))
  
  # Fill matrix row-wise
  for (v in vics) {
    lrs = LR.list[[v]]
    lrs = lrs[names(lrs) != "*"]  # remove do-nothing move
    LR.table[v, names(lrs)] = unname(lrs)
  }
  
  # Reduce moves according to `limit` and/or nkeep
  moves.reduced = lapply(LR.list, function(lrs) {
    newmoves = names(lrs)[lrs > limit]
    if(!is.null(nkeep) && length(newmoves) > nkeep)
      length(newmoves) = nkeep
    newmoves
  })
  
  list(moves = moves.reduced, LR.list = LR.list, LR.table = LR.table)
}
   


# Not used 
screen1 = function(pm, am, missing, moves, loglik0, vict = 1, LRlimit = 0.1, 
                   verbose = F, sorter = FALSE, nkeep = NULL){
  ids.pm = unlist(labels(pm))
  marks = 1:nMarkers(pm)
  idFrom = names(moves[vict])  # id of victim to potentially move
  from1 = pm[[idFrom]] # singleton corresponding to idFrom
  move1 = moves[[vict]] # ids of potential moves
  idTo = as.character(move1)
  nm = length(move1)
  
  # Loglik of each victim
  logliks.PM = vapply(pm, loglikTotal, markers = marks, FUN.VALUE = 1)
  names(logliks.PM) = ids.pm
  
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
      ids.remaining = setdiff(ids.pm, idFrom)
      loglik.remaining = sum(logliks.PM[ids.remaining])
      
      # Likelihood of families after move
      am2 = transferMarkers(pm, am, idsFrom = idFrom, 
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
    list(loglik = logl, LR = LR, keep = keep, move1Kept = move1)
}
  
