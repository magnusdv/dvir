#' Joint DVI search
#'
#' Victims are given as a list of singletons and references as list of
#' pedigrees. Based on a specification of moves for each victim, if any, all
#' possible assignments are evaluated and solutions ranked according to the
#' likelihood.
#'
#' @param pm A list of singletons.
#' @param am A list of pedigrees.
#' @param missing Character vector with names of missing persons.
#' @param moves A list of possible assignments for each victim. If NULL, all
#'   sex-matching assignments are considered.
#' @param limit A positive number. Only single-search LR values above this are
#'   considered.
#' @param markers A vector indicating which markers should be included in the
#'   analysis. By default all markers are included.
#' @param disableMutations A logical, or NA (default). The default action is to
#'   disable mutations in all reference families without Mendelian errors.
#' @param undisputed A logical.
#' @param threshold A positive number, passed onto [findUndisputed()].
#' @param numCores Integer. The number of cores used in parallelisation.
#'   Default: 1.
#' @param check A logical, indicating if the input data should be checked for
#'   consistency.
#' @param verbose A logical.
#'
#' @return A data frame. Each row describes an a priori possible move. The log
#'   likelihood, the LR with the null hypothesis of no moves in the numerator,
#'   and the posterior corresponding to a flat prior, i.e., the inverse of the
#'   number of assignments.
#'
#' @seealso [singleSearch()]
#'
#' @examples
#' pm = example2$pm
#' am = example2$am
#' missing = example2$missing
#' 
#' jointDVI(pm, am, missing)
#' \donttest{
#' # Last example in Chapter 4 of 'Mass identifications", Kling et al (2021)
#'   pm = dataCh4$pm
#'   am = dataCh4$am
#'   missing = dataCh4$missing
#'   plotPedList(am, marker = 1:2)
#'   plotPedList(pm, marker = 1:2)
#'   res = jointDVI(pm, am, missing)
#'   res = jointDVI(pm, am, missing, disableMutations = FALSE)
#'   res[c(1,2,30,49),]
#'   }
#' 
#' @importFrom parallel makeCluster stopCluster detectCores parLapply
#'   clusterEvalQ clusterExport
#'
#' @export
jointDVI = function(pm, am, missing, moves = NULL, limit = 0, undisputed = TRUE, markers = NULL,
                    threshold = 1e4, disableMutations = NA, numCores = 1, check = TRUE, verbose = FALSE){
  
  st = Sys.time()
  
  if(is.singleton(pm)) 
    pm = list(pm)
  if(is.ped(am)) 
    am = list(am)
  
  names(pm) = origVics = vics = unlist(labels(pm)) 
  
  if(!is.null(markers)) {
    pm = selectMarkers(pm, markers)
    am = selectMarkers(am, markers)
  }
  
  if(verbose)
    summariseDVI(pm, am, missing, printMax = 10)
  
  if(check)
    checkDVI(pm, am, missing, moves = moves)
  
  ### Mutation disabling
  if(any(allowsMutations(am))) {
    
    if(verbose) 
      message("\nMutation modelling:")
    
    if(isTRUE(disableMutations)) {
      if(verbose) message(" Disabling mutations in all families")
      disableFams = seq_along(am)
    }
    else if(identical(disableMutations, NA)) {
      am.nomut = pedprobr::setMutationModel(am, model = NULL)
      badFams = vapply(am.nomut, loglikTotal, FUN.VALUE = 1) == -Inf
      if(verbose) {
        if(any(badFams)) 
          message(" ", sum(badFams), " inconsistent families: ", trunc(which(badFams)))
        message(" ", sum(!badFams), " consistent families. Disabling mutations in these")
      }
      disableFams = which(!badFams)
    }
    else disableFams = NULL
  
    if(length(disableFams)) {
      am[disableFams] = setMutationModel(am[disableFams], model = NULL)
    }
  }
  
  ### Identify and fixate "undisputed" matches
  undisp = list()
  
  if(undisputed) {
    
    if(verbose) {
      message("\nUndisputed matches:")
      message(" Single-search LR threshold = ", threshold)
    }
    
    r = findUndisputed(pm, am, missing, moves = moves, threshold = threshold, 
                       limit = limit, check = FALSE, verbose = verbose)
    
    # List of undisputed, and their LR's
    undisp = r$undisp 
    
    # If all are undisputed, return early
    if(length(undisp) == length(pm)) {
      solution = lapply(undisp, function(v) v$match)
      
      # Hack to get consistent output: Run through jointDVI() with the solution as `moves` 
      res = jointDVI(pm, am, missing, moves = solution, undisputed = FALSE,
                     markers = markers, threshold = NULL, check = FALSE, verbose = FALSE)
      return(res)
    }
    
    # Reduced DVI problem to be used in the joint analysis
    pm = r$pmReduced
    am = r$amReduced
    missing = r$missingReduced
    vics = names(pm)
    
    # Moves: These exclude those with LR = 0!
    moves = r$moves
  }
    
  if(is.null(moves)) {
    moves = singleSearch(pm, am, missing = missing, moves = moves, limit = limit)$moves
  }
 
  # Expand moves to grid
  moveGrid = expand.grid.nodup(moves)
  nMoves = nrow(moveGrid)
  if(nMoves == 0)
    stop("No possible solutions!")
  if(verbose)
    message("\nAssignments to consider in the joint analysis: ", nMoves, "\n")
  
  # Convert to list; more handy below
  moveList = lapply(1:nMoves, function(i) as.character(moveGrid[i, ]))
  
  # Initial loglikelihoods
  logliks.PM = vapply(pm, loglikTotal, FUN.VALUE = 1)
  logliks.AM = vapply(am, loglikTotal, FUN.VALUE = 1)
  
  loglik0 = sum(logliks.PM) + sum(logliks.AM)
  if(loglik0 == -Inf)
    stop("Impossible initial data: AM component ", toString(which(logliks.AM == -Inf)))
  
  
  # Parallelise
  if(numCores > 1) {
    cl = makeCluster(numCores)
    on.exit(stopCluster(cl))
    clusterEvalQ(cl, library(dvir))
    clusterExport(cl, "singleMove", envir = environment())
    
    if(verbose) message("Using ", length(cl), " cores")
    
    # Loop through moves
    loglik = parLapply(cl, moveList, function(move) 
      singleMove(pm, am, vics, move, loglik0, logliks.PM, logliks.AM))
  }
  else {
    loglik = lapply(moveList, function(move) 
      singleMove(pm, am, vics, move, loglik0, logliks.PM, logliks.AM))
  }
  
  loglik = unlist(loglik)
  
  LR = exp(loglik - loglik0)
  posterior = LR/sum(LR) # assumes a flat prior
  
  # Add undisputed matches
  if(length(undisp)) {
    # Add ID columns
    for(v in names(undisp)) moveGrid[[v]] = undisp[[v]]$match
    
    # Fix ordering
    moveGrid = moveGrid[origVics]
    
    # Fix LR: Multiply with that of the undisputed
    LR = LR * prod(sapply(undisp, `[[`, "LR"))
  }
    
  # Collect results
  tab = cbind(moveGrid, loglik = loglik, LR = LR, posterior = posterior)
  
  # Remove matches with LR <= limit
  #keep = LR > limit
  #if(length(keep) == 0)
  #  stop("No possible assignments. Try reducing limit")
  #tab = tab[keep, , drop = FALSE]
  
  # Sort in decreasing order
  tab = tab[order(tab$loglik, decreasing = TRUE), , drop = FALSE]
  rownames(tab) = NULL
  
  if(verbose)
    message("Time used: ", format(Sys.time() - st, digits = 3))
  
  tab
}



# Function for computing the total log-likelihood after a given move
singleMove = function(pm, am, vics, move, loglik0, logliks.PM, logliks.AM) {
  
  # Victims which actually move
  vicMove = vics[move != "*"]
  mpsMove = move[move != "*"]
  
  if(length(vicMove) == 0)
    return(loglik0)
  
  # The relevant AM components
  compNo = unique.default(getComponent(am, mpsMove, checkUnique = TRUE))
  
  # Move victim data
  changedComps = transferMarkers(pm[vicMove], am[compNo], idsFrom = vicMove, 
                                 idsTo = mpsMove, erase = FALSE)

  # Update likelihood of modified AM comps
  logliks.AM.new = logliks.AM
  logliks.AM.new[compNo] = vapply(changedComps, function(a) loglikTotal(a), FUN.VALUE = 1)
  
  # Likelihood of remaining PMs
  logliks.PM.new = logliks.PM[setdiff(vics, vicMove)]
  
  # Return total loglik after move
  loglik.move = sum(logliks.PM.new) + sum(logliks.AM.new)
}



# @rdname jointDVI
# @export
checkDVI = function(pm, am, missing, moves){
  if(is.null(moves))
    return()
  # If moves are already expanded, skip checks
  if(is.data.frame(moves))
    return()
  
  if(length(missing) < 1)
    stop("Third vector must a vector with names of missing persons")
  if(!is.list(moves))
    stop("Fourth argument must be a list")
  
  # Check that all specified moves are sex consistent
  idsMoves = names(moves)
  if(any(duplicated(idsMoves)))
    stop("Duplicated names of moves")
  
  sexMoves = getSex(pm, idsMoves)
  for (i in 1:length(moves)){
    candidates = setdiff(moves[[i]], "*")
    if(length(candidates) > 0 ){
      if(!all(candidates %in% missing))
        stop("Wrong mp in element ", i, " of moves")
      thisCheck = all(sexMoves[i] == getSex(am, candidates))
      if(!thisCheck)
        stop("Wrong sex in element ", i, " of moves")
    }
    if(length(candidates) == 0 & idsMoves[[i]] != idsMoves[i])
      stop("Wrong identity move in element ", i, " of moves")
  }
}

  
summariseDVI = function(pm, am, missing, printMax = 10) {
  vics = unlist(labels(pm))
  refs = typedMembers(am)
  
  message("DVI problem:")
  message(sprintf(" %d victims: %s", length(pm), trunc(vics, printMax)))
  message(sprintf(" %d missing: %s", length(missing), trunc(missing, printMax)))
  message(sprintf(" %d families", length(am)))
  message(sprintf(" %d typed refs: %s", length(refs), trunc(refs, printMax)))
}


trunc = function(x, printMax = 10) {
  if(length(x) <= printMax)
    return(toString(x))
  y = c(x[1:5], "...", x[length(x)])
  toString(y)
}
