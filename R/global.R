#' Search for best DVI solution and check input
#' 
#' Victims are given as a list of singletons and references as list of pedigrees.
#' Based on a specification of moves for each victim, if any, all possible assignments 
#' are evaluated and solutions ranked according to the likelihood.. 
#' @aliases global
#' @aliases checkDVI
#' @param from A list of singletons. 
#' @param to A list of pedigrees.
#' @param ids.to Character vector with names of missing persons.
#' @param moves List of length equal length of `from` with possible marginal moves.
#' @param limit Double. Only moves with LR above limit are kept.
#' @param numCores Integer. The number of cores used in parallelisation. Default: 1.
#' @param check A logical, indicating if the input data should be checked for consistency.
#' @param verbose A logical.
#' @details This is currently a brute force approach, all possibilities are evaluated
#' 
#' @return A data frame. Each row describes an a priori possible move. The log likelihood, 
#' the LR with the null hypothesis of no moves in the numerator, and the posterior corresponding 
#' to a flat prior, i.e., the inverse of the number of assignments 
#' @seealso `marginal`.

#' @import pedtools
#' @import forrel
#' @importFrom pedprobr likelihood
#'
#' @examples
#' \donttest{
#' library(forrel)
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
#' from = as.ped(df, locusAttributes = locAttr)
#'
#' # AM data (families)
#' MPs = c("MP1", "MP2", "MP3")
#' to = nuclearPed(3, father = "R1", mother = "R2", children = MPs)
#' m = marker(to, "R1" = "1/1", "R2" = "1/1", name = "m")
#' to = setMarkers(to, m, locusAttributes = locAttr)
#'
#' # Plot
#' plotPedList(list(from, to), marker = 1, col = list(red = MPs),
#'             titles = c("PM data", "AM data"))
#'
#' # Analysis considering all sex-consistent assignments
#' res1 = global(from, to, MPs, moves = NULL, limit = -1, verbose = FALSE)
#'
#' # Quicker alternative: Consider only the three best moves for each victim
#' moves = generateMoves(from, to, MPs) # generate all sex-consistent assignments
#' moves2 = marginal(from, to,  MPs, moves, limit = -1, nkeep = 3)
#' res2 = global(from, to, MPs, moves2[[1]], limit = -1, verbose = FALSE)
#'
#' # Further reduction: Only consider victims V1, V3 and V4
#' moves2 = moves[c("V1", "V3", "V4")]
#' res = global(from, to, MPs, moves2, limit = -1)
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
#' from = as.ped(df, locusAttributes = loc)
#'
#' # Reference families: 4 missing persons; 2 genotyped relatives
#' fam1 = nuclearPed(1, father = "MP1", child = "MP2", mother ="R1")
#' fam2 = halfSibPed(sex1 = 1, sex2 = 2)
#' fam2 = relabel(fam2, c("MO2", "R2", "MO3", "MP3", "MP4"))
#' data = data.frame(m1 = c(R1 = "2/2", R2 = "3/3"))
#' to = setMarkers(list(fam1, fam2), alleleMatrix = data, locusAttributes = loc)
#'
#' # Generate sex-consistent moves
#' ids.to = c("MP1", "MP2", "MP3", "MP4") # user specified
#' moves = generateMoves(from, to, ids.to)
#'
#' # Rank according to likelihood
#' res1 = global(from, to, ids.to, moves = moves, limit = 0, verbose = TRUE)
#' resmarg = marginal(from, to, ids.to, moves = moves, limit = 0,
#'              verbose = TRUE, nkeep = 2)
#' res2 = global(from, to, ids.to, moves = resmarg[[1]], limit = 0, verbose = TRUE)
#'
#' # From fam file
#' x = readFam("http://familias.name/BookKEP/globalExample2.fam")
#' from = x$families[1:4]
#' to = x$families[5:6]
#' plotPedList(list(from, to))
#' res3 = global(from, to, ids.to, moves = moves, limit = 0, verbose = TRUE)
#'
#' stopifnot(identical(res1, res3))
#'
#' }
#'
#' @importFrom parallel makeCluster stopCluster detectCores parLapply
#'   clusterEvalQ
#' @export global checkDVI



global = function(from, to, ids.to, moves = NULL, limit = 0, numCores = 1, check = TRUE, verbose = FALSE){
  if(is.null(moves)) # Generate assignments
    moves = generateMoves(from = from, to = to, ids.to = ids.to, expand.grid = TRUE)
  
  # Check consistency
  if(check)
    checkDVI(from = from, to = to, ids.to = ids.to, moves = moves)
  
  ids.from = unlist(labels(from))  # as.character(lapply(from, function(x) x$ID))
  marks = seq_len(nMarkers(from))  # 1:nMarkers(from)
  
  # Initial loglikelihoods
  logliks.PM = vapply(from, loglikTotal, markers = marks, FUN.VALUE = 1)
  names(logliks.PM) = ids.from
  
  loglik0 = sum(logliks.PM) + loglikTotal(to, markers = marks)
  if(loglik0 == -Inf)
    stop("Impossible initial data")
  
  # Expand if needed
  moveGrid = if(is.data.frame(moves)) moves else expand.grid.nodup(moves)
  nMoves = nrow(moveGrid)
  if(nMoves == 0)
    stop("No possible solutions specified, possibly identical moves")

  vics = names(moveGrid)
  
  # Convert to list, which is more handy below
  moveList = lapply(1:nMoves, function(i) as.character(moveGrid[i, ]))
  
  # Function for computing the total log-likelihood after a given move
  singleMove = function(from, to, vics, move, loglik0, logliks.PM) {
    
    # Victims which actually move
    move.from = vics[move != "*"]
    move.to = move[move != "*"]
    
    if(length(move.from) == 0)
      return(loglik0)
   
    # Likelihood of remaining PMs
    remaining = setdiff(vics, move.from)
    loglik.remaining = sum(logliks.PM[remaining])
    
    # Likelihood of families after move
    to2 = transferMarkers(from, to, idsFrom = move.from,
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
      singleMove(from, to, vics, move, loglik0, logliks.PM))
  }
  else {
    loglik = lapply(moveList, function(move) 
      singleMove(from, to, vics, move, loglik0, logliks.PM))
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
  
  tab
}

checkDVI = function(from, to, ids.to, moves){
  # If moves are already expanded, skip checks
  if(is.data.frame(moves))
    return()
  
  if(is.null(moves))
    stop("No moves specified")
  if(length(ids.to) < 1)
    stop("Third vector must a vector with names of missing persons")
  if(!is.list(moves))
    stop("Fourth argument must be a list")
  
  # Check that all specified moves are sex consistent
  idsMoves = names(moves)
  if(any(duplicated(idsMoves)))
    stop("Duplicated names of moves")
  
  sexMoves = getSex(from, idsMoves)
  for (i in 1:length(moves)){
    candidates = setdiff(moves[[i]], "*")
    if(length(candidates) > 0 ){
      if(!all(candidates %in% ids.to))
        stop("Wrong mp in element ", i, " of moves")
      thisCheck = all(sexMoves[i] == getSex(to, candidates))
      if(!thisCheck)
        stop("Wrong sex in element ", i, " of moves")
    }
    if(length(candidates) == 0 & idsMoves[[i]] != idsMoves[i])
      stop("Wrong identity move in element ", i, " of moves")
  }
  
}

  
