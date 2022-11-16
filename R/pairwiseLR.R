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
#' @param limit A nonnegative number controlling the `pairing` slot of the
#'   output: Only pairings with LR greater or equal to `limit` are kept. If zero
#'   (default), pairings with LR > 0 are kept.
#' @param nkeep An integer. No of pairings to keep, all if `NULL`.
#' @param check A logical, indicating if the input data should be checked for
#'   consistency.
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
#'
#' @examples
#' pairwiseLR(example1, verbose = TRUE)
#'
#' @export
pairwiseLR = function(dvi, pairings = NULL, ignoreSex = FALSE, limit = 0, nkeep = NULL, 
                    check = TRUE, verbose = FALSE){
  
  if(verbose)
    message("Computing matrix of pairwise LR")
  
  # Ensure proper dviData object
  dvi = consolidateDVI(dvi)
  
  # Check consistency
  if(check)
    checkDVI(dvi, pairings = pairings, ignoreSex = ignoreSex)
  
  pm = dvi$pm
  am = dvi$am
  missing = dvi$missing
  
  if(length(pm) == 0)
    return(list(LRmatrix = NULL, LRlist = list(), pairings = list()))
  
  if(is.null(pairings)) # Generate pairings
    pairings = generatePairings(dvi, ignoreSex = ignoreSex)
  
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
  LRlist = lapply(vics, function(v) {
    
    # Corresponding vector of LRs
    lrs = vapply(pairings[[v]], function(mp) {
      
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
      
      # Return LR
      exp(loglik.move - loglik0)
    }, FUN.VALUE = numeric(1))
    
    # Return sorted vector
    sort(lrs, decreasing = TRUE)
  })
  
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
  pairings.reduced = lapply(LRlist, function(lrs) {
    newpairings = names(lrs)[lrs > 0 & lrs >= limit]
    if(!is.null(nkeep) && length(newpairings) > nkeep)
      length(newpairings) = nkeep
    newpairings
  })
  
  list(LRmatrix = LRmatrix, LRlist = LRlist, pairings = pairings.reduced)
}

