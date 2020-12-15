#' Check and evaluate individual candidates
#' 
#' Each potential move is evaluated and the search space reduced. 
#' 
#' @param pm A list of singletons, the victims. 
#' @param am A list of pedigrees. The reference families.
#' @param missing Character vector with names of missing persons.
#' @param moves List with possible marginal moves.
#' @param limit Double. Lower threshold for LR.
#' @param nkeep integer. No of moves to keep, all if `NULL`.
#' @param check A logical, indicating if the input data should be checked for consistency.
#' @param verbose Logical.
#' @details The potential reduction only affects the list of moves returned, all LRs are kept.
#' Specifying `nkeep`can give further reduction.
#' 
#' 
#' @return A list with moves and log likelihoods.
#' 
#' @examples
#' 
#' \donttest{
#' 
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
#' pm = as.ped(df, locusAttributes = loc)
#' names(pm) = vics
#' 
#' # Reference families
#' missing = c("MP1", "MP2", "MP3")
#' am = nuclearPed(3, father = "R1", mother = "R2", children = missing)
#' data = data.frame(m = c("1/1", "1/1"), row.names = c("R1", "R2"))
#' am = setMarkers(am, alleleMatrix = data, locusAttributes = loc)
#' 
#' plotPedList(list(pm, am), marker = 1)
#' 
#' # moves = list(V1 = c("MP1", "V1", "MP2"), 
#' #              V3 = c("MP1", "MP2", "MP3"), 
#' #              V4 = c("V4", "MP3"), 
#' #              V7 = c("V7"))
#' #             
#' # # all 4 * 6 + 1 =possible moves ignoring sex
#' # moves = list(V1 = c("V1", missing), V2 = c("V2", missing), V3 = c("V3", missing),
#' #              V4 = c("V4", missing), V5 = c("V5", missing), V6 = c("V6", missing),
#' #              V7 = c("V7"))
#' #
#' # res = marginal(pm, am, missing, moves, limit = 0, verbose = TRUE, nkeep= 2)
#'  
#' moves = generateMoves(pm, am, missing)
#' 
#' res = marginal(pm, am, missing, moves, limit = -1, nkeep = 3)
#' res2 = jointDVI(pm, am, missing, moves = res[[1]], limit = 0)
#' # moves = list(V1 = c("V1", "MP1", "MP2"))
#' # res = marginal(pm, am,  missing, moves, limit = 1)
#' }
#' 
#' 
#' @export
marginal = function(pm, am, missing, moves = NULL, limit = 0.1, nkeep = NULL, 
                    check = TRUE, verbose = FALSE){
  
  if(is.singleton(pm)) pm = list(pm)
  if(is.ped(am)) am = list(am)
  
  if(is.null(moves)) # Generate moves
    moves = generateMoves(pm = pm, am = am, missing = missing)
  
  # Check consistency
  if(check)
    checkDVI(pm = pm, am = am, missing = missing, moves = moves)

  marks = 1:nMarkers(pm)
  
  # Ensure correct names
  names(pm) = unlist(labels(pm), use.names = FALSE)
  vics = names(moves) # normally in the same order as names(pm)
  
  # Loglik of each victim
  logliks.PM = vapply(pm, loglikTotal, markers = marks, FUN.VALUE = 1)
  
  # Loglik of each ref family
  logliks.AM = vapply(am, loglikTotal, markers = marks, FUN.VALUE = 1)
  
  # log-likelihood of H0
  loglik0 = sum(logliks.PM) + sum(logliks.AM)
  
  if(loglik0 == -Inf)
    stop("Impossible initial data")
  
  # For each victim, compute marginal LRs of each move
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
  
  # Matrix of marginal LRs (filled with 0's)
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
    newmoves = names(lrs)[lrs >= limit]
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
    list(loglik = logl, LR = LR, keep = keep, 
              move1Kept = move1)
}
  


