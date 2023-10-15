#' Pairwise LR matrix
#'
#' For a given DVI problem, compute the matrix consisting of pairwise likelihood
#' ratios \eqn{LR_{i,j}} comparing \eqn{V_i = M_j} to the null. The output may
#' be reduced by specifying arguments `limit` or `nkeep`.
#'
#' @param dvi A `dviData` object, typically created with [dviData()].
#' @param pairings A list of possible pairings for each victim. If NULL, all
#'   sex-consistent pairings are used.
#' @param ignoreSex A logical.
#' @param limit A nonnegative number controlling the `pairings` slot of the
#'   output: Only pairings with LR greater or equal to `limit` are kept. If zero
#'   (default), pairings with LR > 0 are kept.
#' @param nkeep An integer, or NULL. If given, only the `nkeep` most likely
#'   pairings are kept for each victim.
#' @param check A logical, indicating if the input data should be checked for
#'   consistency.
#' @param numCores An integer; the number of cores used in parallelisation.
#'   Default: 1.
#' @param verbose A logical.
#'
#' @return A list with 3 elements:
#'
#'   * `LRmatrix`: A matrix containing the pairwise LR values.
#'
#'   * `LRlist`: A list of numerical vectors, containing the pairwise LRs in
#'   list format.
#'
#'   * `pairings`: A reduced version of the input `pairings`, keeping only
#'   entries with corresponding LR >= `limit`. For the default case `limit = 0`
#'   a strict inequality is used, i.e., LR > 0.
#'
#' @examples
#' pairwiseLR(example1, verbose = TRUE)
#'
#' @export
pairwiseLR = function(dvi, pairings = NULL, ignoreSex = FALSE, limit = 0, nkeep = NULL, 
                    check = TRUE, numCores = 1, verbose = FALSE){
  
  if(verbose)
    cat("Computing matrix of pairwise LR\n")
  
  # Ensure proper dviData object
  dvi = consolidateDVI(dvi)
  
  # Check consistency
  if(check)
    checkDVI(dvi, pairings = pairings, ignoreSex = ignoreSex, verbose = verbose)
  
  pm = dvi$pm
  am = dvi$am
  missing = dvi$missing
  
  if(length(pm) == 0)
    return(list(LRmatrix = NULL, LRlist = list(), pairings = list()))
  
  # Generate pairings
  pairings = pairings %||% dvi$pairings %||% generatePairings(dvi, ignoreSex = ignoreSex)
  
  # Loglik of each victim and each ref
  marks = 1:nMarkers(pm)
  logliks.PM = vapply(pm, loglikTotal, markers = marks, FUN.VALUE = 1)
  logliks.AM = vapply(am, loglikTotal, markers = marks, FUN.VALUE = 1)
  
  # log-likelihood of H0
  loglik0 = sum(logliks.PM) + sum(logliks.AM)
  
  if(loglik0 == -Inf)
    stop2("Impossible initial data: AM component ", which(logliks.AM == -Inf))
  
  # For each victim, compute the LR of each pairing
  vics = names(pm)
  
  # Dont use more cores than the number of vics
  numCores = min(numCores, length(vics))
  
  # Parallelise
  if(numCores > 1) {
    
    if(verbose) 
      cat("Using", numCores, "cores\n")
    
    cl = makeCluster(numCores)
    on.exit(stopCluster(cl))
    clusterEvalQ(cl, library(dvir))
    clusterExport(cl, "pairwise_singlevic", envir = environment())
    
    # Loop through victims
    LRlist = parLapply(cl, vics, function(v) 
      pairwise_singlevic(am, vics, v, pm[[v]], pairings[[v]], marks, loglik0, logliks.PM, logliks.AM))
  }
  else {
    # Default: no parallelisation
    LRlist = lapply(vics, function(v) 
      pairwise_singlevic(am, vics, v, pm[[v]], pairings[[v]], marks, loglik0, logliks.PM, logliks.AM))
  }
  
  names(LRlist) = vics
  
  # Matrix of individual LRs (filled with 0's)
  LRmatrix = matrix(0, nrow = length(vics), ncol = length(missing), 
                    dimnames = list(vics, missing))
  
  # Fill matrix row-wise
  for (v in vics) {
    lrs = LRlist[[v]]
    lrs = lrs[names(lrs) != "*"]  # remove do-nothing move
    LRmatrix[v, names(lrs)] = unname(lrs)
  }
  
  # Reduce pairings according to `limit` and/or nkeep
  pairingsReduced = lapply(LRlist, function(lrs) {
    if(limit == 0) 
      keepIdx = lrs > 0
    else
      keepIdx = lrs >= limit
    if(limit > 1) 
      keepIdx = keepIdx | names(lrs) == "*"  # Never remove "*"
    
    newpairings = names(lrs)[keepIdx]
    
    # Apply `nkeep` if given
    if(!is.null(nkeep) && length(newpairings) > nkeep) {
      starIdx = match("*", newpairings, nomatch = 0)
      length(newpairings) = nkeep
      if(starIdx > nkeep)
        newpairings = c(newpairings, "*")
    }
    
    newpairings
  })
  
  list(LRmatrix = LRmatrix, LRlist = LRlist, pairings = pairingsReduced)
}


# Function for computing the pairwise LRs for a single victim
pairwise_singlevic = function(am, vics, v, pmV, pairingsV, marks, loglik0, logliks.PM, logliks.AM) {
  
  lrs = vapply(pairingsV, function(mp) {
    
    if(mp == "*") 
      return(1)
    
    # Make copy of AM likelihoods (vector)
    logliks.AM.new = logliks.AM
    
    # The relevant AM component 
    compNo = getComponent(am, mp, checkUnique = TRUE)
    
    # Move victim data to `mp`
    comp = transferMarkers(pmV, am[[compNo]], idsFrom = v, idsTo = mp, erase = FALSE)
    
    # Update likelihood of this comp
    logliks.AM.new[compNo] = loglikTotal(comp, marks)
    
    # Likelihood of remaining PMs
    logliks.PM.new = logliks.PM[setdiff(vics, v)]
    
    # Total loglik after move
    loglik.move = sum(logliks.PM.new) + sum(logliks.AM.new)
    
    # Return LR
    exp(loglik.move - loglik0)
  }, FUN.VALUE = numeric(1))
  
  # Return sorted vector
  sort(lrs, decreasing = TRUE)
}