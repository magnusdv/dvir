#' Joint DVI search
#' 
#' Victims are given as a list of singletons and references as list of pedigrees.
#' Based on a specification of moves for each victim, if any, all possible assignments 
#' are evaluated and solutions ranked according to the likelihood.. 
#' @aliases global
#' @aliases checkDVI
#' @param pm A list of singletons. 
#' @param am A list of pedigrees.
#' @param missing Character vector with names of missing persons.
#' @param moves List of length equal length of `pm` with possible marginal moves.
#' @param limit Double. Only moves with LR above limit are kept.
#' @param fixUndisputed A logical.
#' @param threshold A positive number, passed onto [findUndisputed()].
#' @param numCores Integer. The number of cores used in parallelisation. Default: 1.
#' @param check A logical, indicating if the input data should be checked for consistency.
#' @param verbose A logical.
#' 
#' @return A data frame. Each row describes an a priori possible move. The log likelihood, 
#' the LR with the null hypothesis of no moves in the numerator, and the posterior corresponding 
#' to a flat prior, i.e., the inverse of the number of assignments.
#' 
#' @seealso `marginal`.
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
#' missing = c("MP1", "MP2", "MP3")
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
#' moves2 = marginal(pm, am,  missing, moves, limit = -1, nkeep = 3)
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
#' missing = c("MP1", "MP2", "MP3", "MP4") # user specified
#' moves = generateMoves(pm, am, missing)
#'
#' # Rank according to likelihood
#' res1 = jointDVI(pm, am, missing, moves = moves, limit = 0, verbose = TRUE)
#' resmarg = marginal(pm, am, missing, moves = moves, limit = 0,
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
#'   clusterEvalQ
#'
#' @export
jointDVI = function(pm, am, missing, moves = NULL, limit = 0, fixUndisputed = TRUE, 
                    threshold = 1e4, numCores = 1, check = TRUE, verbose = FALSE){
  
  st = Sys.time()
  
  if(is.singleton(pm))
    pm = list(pm)
  
  vics = unlist(labels(pm)) 
  names(pm) = vics
  
  if(verbose) {
    message("Input data:")
    message(" Victims: ", toString(vics))
    message(" Missing: ", toString(missing))
    message(" Families: ", if(is.ped(am)) 1 else length(am))
    message(" Refs: ", toString(typedMembers(am)))
  }
  
  if(check)
    checkDVI(pm, am, missing, moves = moves)
  
  if(fixUndisputed) {
    undisp = findUndisputed(pm, am, missing, threshold = threshold, check = FALSE, verbose = verbose)
    missing = setdiff(missing, undisp)
    
    # Remaining pm data
    pmRem = pm[setdiff(vics, names(undisp))]
  
    # Generate all moves with the remaining victims
    movesRem = generateMoves(pmRem, am, missing, expand.grid = FALSE)
    
    # Add the undisputed and reorder
    moves = c(movesRem, as.list(undisp))[vics]
  }
    
  if(is.null(moves)) # Generate assignments
    moves = generateMoves(pm = pm, am = am, missing = missing, expand.grid = TRUE)
  
  # Expand if needed
  moveGrid = if(is.data.frame(moves)) moves else expand.grid.nodup(moves)
  nMoves = nrow(moveGrid)
  if(nMoves == 0)
    stop("No possible solutions specified, possibly identical moves")
  
  # Convert to list, which is more handy below
  moveList = lapply(1:nMoves, function(i) as.character(moveGrid[i, ]))
  
  if(verbose) {
    if(fixUndisputed) {
      if(length(undisp) == 0) message("\nUndisputed: None\n")
      else message("\nUndisputed:\n", sprintf(" %s = %s\n", names(undisp), undisp))
    }
    
    message("Assignments to consider: ", nMoves)
    message("")
  }
  
  marks = seq_len(nMarkers(pm))
  
  # Initial loglikelihoods
  logliks.PM = vapply(pm, loglikTotal, markers = marks, FUN.VALUE = 1)
  names(logliks.PM) = vics
  
  loglik0 = sum(logliks.PM) + loglikTotal(am, markers = marks)
  if(loglik0 == -Inf)
    stop("Impossible initial data")
  
  
  # Function for computing the total log-likelihood after a given move
  singleMove = function(pm, am, vics, move, loglik0, logliks.PM) {
    
    # Victims which actually move
    move.from = vics[move != "*"]
    move.to = move[move != "*"]
    
    if(length(move.from) == 0)
      return(loglik0)
   
    # Likelihood of remaining PMs
    remaining = setdiff(vics, move.from)
    loglik.remaining = sum(logliks.PM[remaining])
    
    # Likelihood of families after move
    to2 = transferMarkers(pm, am, idsFrom = move.from,
                          idsTo = move.to, erase  = FALSE)
    loglik.fam = loglikTotal(to2)
    
    # Return total
    loglik.remaining + loglik.fam
  }
  
  
  # Parallelise
  if(numCores > 1) {
    cl = makeCluster(numCores)
    on.exit(stopCluster(cl))
    clusterEvalQ(cl, library(dvir))
    
    if(verbose) message("Using ", length(cl), " cores")
    
    # Loop through moves
    loglik = parLapply(cl, moveList, function(move) 
      singleMove(pm, am, vics, move, loglik0, logliks.PM))
  }
  else {
    loglik = lapply(moveList, function(move) 
      singleMove(pm, am, vics, move, loglik0, logliks.PM))
  }
  
  loglik = unlist(loglik)
  
  LR = exp(loglik - loglik0)
  posterior = LR/sum(LR) # assumes a flat prior
  
  # Collect results
  tab = cbind(moveGrid, loglik = loglik, LR = LR, posterior = posterior)
  
  # Remove matches with LR <= limit
  keep = LR > limit
  if(length(keep) == 0)
    stop("No possible assignments. Try reducing limit")
  tab = tab[keep, , drop = FALSE]
  
  # Sort in decreasing order
  tab = tab[order(tab$loglik, decreasing = TRUE), , drop = FALSE]
  rownames(tab) = NULL
  
  if(verbose)
    message("Time used: ", format(Sys.time() - st, digits = 3))
  
  tab
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

  
