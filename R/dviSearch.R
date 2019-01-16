#' Search for DVI solution
#' 
#' A solution is obtain by finding the most likely among marginal 
#' possible moves
#' 
#' @param from A list of singeltons 
#' @param to A list of pedigrees
#' @param ids.from Character vector with names of victims
#' @param ids.to Character vector with names of missing persons
#' @param limit Double. Threshold for LR
#' @param nbest Integer If not `NULL`, search is restricted to the `nbest`
#' marginal solutions for each victim
#' @param extend Logical. If `TRUE`any number of victims need not be among 
#' the missing persons
#' @param pl Logical. If `TRUE`, a plot of the best solution will be made. A plot window should 
#' be prepared and sized in advance
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
#'      nbest = NULL, extend = T)
#' date()
#' }

dviSearch = function(from, to, ids.from, ids.to, limit = 1, 
                     nbest = 2, extend = TRUE, pl = FALSE){
  f = function(x){
    rownames(g) = x
    am2 = setAlleles(to2, ids = x, alleles = g)
    pm2 = setAlleles(from, ids = ids.from, alleles = 0)
    prod(LR(list(pm2, am2), 1)$likelihoodsPerSystem)
  }
    
  #Check and and find null likelihood
  if(!is.null(nbest))
    if(nbest < 1)
      stop("nbest must be NULL or an integer greater than 0")
  check = checkInput(from, to, ids.from, ids.to)
  if(!is.null(check$error)) stop(check$error)
  lik0 = check$lik0
  #Evaluate marginal moves
  moves = reduce(from, to, ids.from, ids.to, limit = limit)
  if(is.null(moves)){
    sortedMoves = data.frame(matrix(ids.from, nrow = 1), lik = lik0, 
                             lik0 = lik0, LR = 1, row.names = 1)
    colnames(sortedMoves) = c(ids.from, "lik", "lik0", "LR")
    if(pl)
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
  liks = unlist(lapply(moves, f))
  # liks = rep(NA, n.moves)
  # for (i in 1:n.moves){
  #   rownames(g) = moves[[i]]
  #   am2 = setAlleles(to2, ids = moves[[i]], alleles = g)
  #   pm2 = setAlleles(from, ids = ids.from, alleles = 0)
  #   liks[i] = prod(LR(list(pm2, am2), 1)$likelihoodsPerSystem)
  # }
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
  if(pl){   # plot for best solution
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
