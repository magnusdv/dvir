#' Search for DVI solution
#' This is the main function for Disaster Victim Identification (DVI) application
#' implemented in the [dvir] library.
#' 
#' The data consists of post mortem (PM) and ante mortem (AM) samples.
#' The PM samples may be degraded and drop-outs could occur, and possibly also
#' other artefact like drop-in and genotyping error. 
#' There may be samples for one individal and such samples should be merged.
#' Merging can be done in several ways. For instance a consensus profile 
#' amounts to keeping marker data that coincides for all markers. In addition
#' to replicates from one individual, there could also be family relations.
#' The AM data consists of unrelated reference families. Each family has at least on missing person
#' (MP) and at least on genotyped family member. Obviously, the MP-s are not genotyped.
#' We assume that the marker data is of good quality and it is therefore not account for
#' the artefacts mentioned above for the PM data. It is important to model mutations as
#' there will be many comparisons between PM and AM samples.
#' 
#' The main steps of the code are
#' 1. Check input (obvious checking omitted). There should be marker data for all
#' PM samples and no data for MP-s. Each reference family needs at least on genotyped 
#' family member. The missing persons need to be provided as input as there may be other
#' untyped family members included to define family relationships. If initial, null,
#' likelihood is 0, computations are terminated. This typically does not happen 
#' if mutations are modelled.
#' 
#' 2. [Not yet implemented] Identical PM samples and close relations 
#' are identified if the LR compared to
#' unrelatednes exceeds a specified threshold, default 10000. 
#' 
#' 3. 
#' 
#' A solution is obtain by finding the most likely among marginal 
#' possible moves
#' 
#' @param from A list of singeltons 
#' @param to A list of pedigrees
#' @param ids.from Character vector with names of PM-samples if these
#' come from different individuals. If NULL, different individuals are not assumed
#' and names are taken from `from`
#' @param ids.to Character vector with names of missing persons
#' @param limit Double. Threshold for LR
#' @param nbest Integer If not `NULL`, search is restricted to the `nbest`
#' marginal solutions for each victim
#' @param extend Logical. If `TRUE`any number of victims need not be among 
#' the missing persons
#' @param plots Logical. If `TRUE`, plot(s) of the best solution will be made. 
#' A plot window should  #' be prepared and sized in advance
#' @details The search can be restricted using the `nbest` variable.
#' 
#' @return A data frame with moves sorted according to LR. Possibly a plot of the best solution.
#' @export
#' @examples
#' \dontrun{
#' library(forrel)
#' data(dvi.nfi)
#' data(dvi.nfi.mut)
#' dvi.nfi = dvi.nfi.mut
#' from = dvi.nfi$pm
#' to = dvi.nfi$am
#' ids.from = dvi.nfi$vict
#' ids.to = dvi.nfi$miss
#' date()
#' foo = dviSearch(from, to, ids.from,ids.to, limit = -1, 
#'      nbest = 1, extend = T)
#' date()
#' }

dviSearch = function(from, to, ids.from, ids.to, limit = 1, 
                     nbest = 2, extend = TRUE, plots = FALSE){

  # If ids.from is specied, different individuals are assumed, 
  # otherwise sample names are found
  if(is.null(ids.from)) 
    ids.from = getIDs(from)
  
  # Input is checked and NULL likelihood found
  check = checkInput(from, to, ids.from, ids.to, limit = limit, 
                     nbest = nbest, extend = extend, plots = plots)
  if(!is.null(check$error)) stop(check$error)
  lik0 = check$lik0

# Evaluate marginal moves
  moves = reduce(from, to, ids.from, ids.to, limit = limit)
  if(is.null(moves)){
    sortedMoves = data.frame(matrix(ids.from, nrow = 1), lik = lik0, 
                             lik0 = lik0, LR = 1, row.names = 1)
    colnames(sortedMoves) = c(ids.from, "lik", "lik0", "LR")
    if(plots)
      plotPedList(list(from, to), newdev = F, marker = 1,
                  skip.empty.genotypes = TRUE)
    return(sortedMoves)
  }
  if(dim(moves)[1] == 1)
    return(moves)
  moves = moves[,1:2]
  moves = as.data.frame(apply(moves,2, function(x) as.character(x)))
  ids.from = as.character(unique(moves[,1]))
  ids.to = as.character(unique(moves[,2]))
  if(!extend & length(ids.from) > length(ids.to))
    stop("No solution. Try `extend `= TRUE")
  
  #Rearrange input table of moves
  moves = split(moves$to, moves$from)
  moves = lapply(moves, as.character)
  if(!is.null(nbest))
    moves = lapply(moves, function(x) x[1:min(nbest, length(x))])
  to2 = to
  if (extend){ #Add victims to allow no moves
    for (i in 1:length(moves))
      moves[[i]] = c(moves[[i]],ids.from[i])
    ids.to = c(ids.to, ids.from)
    from0 = setAlleles(from, alleles = 0)
    to2 = list()
    if(is.ped(to)){
      n = 1
      to2[[1]] = to
    }
    if(is.pedList(to)){
      n = length(to)
      for(i in 1:n)
        to2[[i]] = to[[i]]
    }
    if(is.ped(from0))
      to2[[n+1]] = from0
    if(is.pedList(from0)){
      n2 = length(from0)
      for(i in (n+1):(n2+n))
        to2[[i]] = from0[[i-n]]
    }
  }
  # Find all paths and corresponding likelihoods
  moves = expand.grid.nodup(moves)
  n.moves = length(moves)
  g = getAlleles(from, ids.from)
  liks = rep(NA, n.moves)
  for (i in 1:n.moves){
    rownames(g) = moves[[i]]
    am2 = setAlleles(to2, ids = moves[[i]], alleles = g)
    pm2 = setAlleles(from, ids = ids.from, alleles = 0)
    liks[i] = prod(LR(list(pm2, am2), 1)$likelihoodsPerSystem)
  }
  lr = liks/lik0
  index = (1:n.moves)[lr > limit]
  lr = lr[index]
  liks = liks[index]
  movesMatrix = matrix(unlist(moves[index]), ncol = length(ids.from), byrow = TRUE)
  ord = order(lr, decreasing = TRUE)
  movesFrame= as.data.frame(movesMatrix)
  sortedMoves = data.frame(movesFrame[ord,], lik = liks[ord], lik0 = lik0, LR = lr[ord],
                             row.names = 1:length(index))
  colnames(sortedMoves) = c(ids.from, "lik", "lik0", "LR")
  if(plots){   # plot for best solution
    best = as.matrix(sortedMoves[1, 1:length(ids.from)], ncol = 1)
    imove = (1:dim(best)[2])[!colnames(best)%in%best]
    g = getAlleles(from, ids.from[imove])
    from = setAlleles(from, ids = ids.from[imove], alleles = 0)
    rownames(g) = best[1,imove]
    to = setAlleles(to, ids = rownames(g), alleles = g)
    plotPedList(list(from, to), newdev = F, marker = 1,
      skip.empty.genotypes = TRUE)
  }
  sortedMoves
}
