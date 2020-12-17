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
#' @param moves List of length equal length of `pm` with possible individual
#'   moves.
#' @param limit A positive number. Only single-search LR values above this are
#'   considered.
#' @param markers A vector indicating which markers should be included in the
#'   analysis. By default all markers are included.
#' @param disableMutations A logical, or NA (default). The default action is to
#'   disable mutations in all reference families without Mendelian errors.
#' @param fixUndisputed A logical.
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
#'
#' \donttest{
#'
#' library(pedtools)
#'
#' ### Example 1 ###
#'
#' # Attributes of a single marker
#' locAttr = list(name = "m", alleles = 1:3, afreq = c(1, 1, 1)/3)
#'
#' # PM data (victims: 6 males, 1 female)
#' vics = paste0("V", 1:7)
#' sex = c(1,1,1,1,1,1,2)
#' df = data.frame(famid = vics, id = vics, fid = 0, mid = 0, sex = sex,
#'                 m = c("1/1", "2/2", "1/1", "1/1", "2/2", "2/2", "2/2"))
#' pm = as.ped(df, locusAttributes = locAttr)
#'
#' # AM data (families)
#' missing = c("M1", "M2", "M3")
#' am = nuclearPed(3, father = "R1", mother = "R2", children = missing)
#' m = marker(am, "R1" = "1/1", "R2" = "1/1", name = "m")
#' am = setMarkers(am, m, locusAttributes = locAttr)
#'
#' # Plot
#' plotPedList(list(pm, am), marker = 1, col = list(red = missing),
#'             titles = c("PM data", "AM data"))
#'
#' # Analysis considering all sex-consistent assignments
#' res1 = jointDVI(pm, am, missing, moves = NULL, limit = -1, verbose = FALSE)
#'
#' # Quicker alternative: Consider only the three best moves for each victim
#' moves = generateMoves(pm, am, missing) # generate all sex-consistent assignments
#' moves2 = singleSearch(pm, am,  missing, moves, limit = -1, nkeep = 3)
#' res2 = jointDVI(pm, am, missing, moves2[[1]], limit = -1, verbose = FALSE)
#'
#' # Further reduction: Only consider victims V1, V3 and V4
#' # moves2 = moves[c("V1", "V3", "V4")]
#' # res = jointDVI(pm, am, missing, moves2, limit = -1)
#'
#' ### Example 2 ###
#' # Checked against; http://familias.name/BookKEP/globalExample2.fam
#'
#' # Single locus with 10 alleles
#' loc = list(alleles = 1:10)
#'
#' # Victims: 3 males, 1 female
#' vics = paste0("V", 1:4)
#' sex = c(1, 1, 2, 1)
#' df = data.frame(famid = vics, id = vics, fid = 0, mid = 0, sex = sex,
#'                 m1 = c("1/1", "1/2", "3/4", "3/4"))
#' pm = as.ped(df, locusAttributes = loc)
#'
#' # Reference families: 4 missing persons; 2 genotyped relatives
#' fam1 = nuclearPed(1, father = "MP1", child = "MP2", mother ="R1")
#' fam2 = halfSibPed(sex1 = 1, sex2 = 2)
#' fam2 = relabel(fam2, c("MO2", "R2", "MO3", "MP3", "MP4"))
#' data = data.frame(m1 = c(R1 = "2/2", R2 = "3/3"))
#' am = setMarkers(list(fam1, fam2), alleleMatrix = data, locusAttributes = loc)
#'
#' # Generate sex-consistent moves
#' missing = paste0("MP", 1:4)
#' moves = generateMoves(pm, am, missing)
#'
#' # Rank according to likelihood
#' res1 = jointDVI(pm, am, missing, moves = moves, limit = 0, verbose = TRUE)
#' resmarg = singleSearch(pm, am, missing, moves = moves, limit = 0,
#'              verbose = TRUE, nkeep = 2)
#' res2 = jointDVI(pm, am, missing, moves = resmarg[[1]], limit = 0, verbose = TRUE)
#'
#' # Victims fam file
#' x = forrel::readFam("http://familias.name/BookKEP/globalExample2.fam")
#' pm = x$families[1:4]
#' am = x$families[5:6]
#' plotPedList(list(pm, am))
#' res3 = jointDVI(pm, am, missing, moves = moves, limit = 0, verbose = TRUE)
#'
#' stopifnot(identical(res1, res3))
#'
#' }
#'
#' @importFrom parallel makeCluster stopCluster detectCores parLapply
#'   clusterEvalQ clusterExport
#'
#' @export
jointDVI = function(pm, am, missing, moves = NULL, limit = 0, fixUndisputed = TRUE, markers = NULL,
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
  
  if(fixUndisputed) {
    
    if(verbose) {
      message("\nUndisputed matches:")
      message(" Single-search LR threshold = ", threshold)
    }
    
    r = findUndisputed(pm, am, missing, threshold = threshold, 
                       limit = limit, check = FALSE, verbose = verbose)
    
    # List of undisputed, and their LR's
    undisp = r$undisp 
    
    # Reduced DVI problem to be used in the joint analysis
    pm = r$pmReduced
    am = r$amReduced
    missing = r$missingReduced
    vics = names(pm)
    
    # Moves: These exclude those with LR = 0!
    moves = r$moves
  }
    
  if(is.null(moves)) {
    moves = singleSearch(pm, am, missing = missing, limit = limit)$moves
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



#' @rdname jointDVI
#' @export
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