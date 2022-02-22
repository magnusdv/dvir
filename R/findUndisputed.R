#' Undisputed identifications in DVI problems
#'
#' This function uses the pairwise LR matrix to find "undisputed" matches
#' between victims and missing individuals. An identification \eqn{V_i = M_j} is
#' called undisputed if the corresponding likelihood ratio \eqn{LR_{i,j}}
#' exceeds the given `threshold`, while all other pairwise LRs involving
#' \eqn{V_i} or \eqn{M_j} are at most 1.
#'
#' If the parameter `relax` is set to TRUE, the last criterion is relaxed,
#' requiring instead that \eqn{LR_{i,j}} is at least `threshold` times greater
#' than all other pairwise LRs involving \eqn{V_i} or \eqn{M_j}
#'
#' @param pm PM data: List of singletons.
#' @param am AM data: A `ped` object or list of such.
#' @param missing Character vector with names of the missing persons.
#' @param pairings A list of possible pairings for each victim. If NULL, all
#'   sex-consistent pairings are used.
#' @param ignoreSex A logical.
#' @param threshold A non-negative number. If no pairwise LR exceed this, the
#'   iteration stops.
#' @param relax A logical affecting the definition of being undisputed (see
#'   Details). Default: FALSE.
#' @param limit A positive number. Only pairwise LR values above this are
#'   considered.
#' @param check A logical indicating if the input data should be checked for
#'   consistency. Default: TRUE.
#' @param verbose A logical. Default: TRUE.
#'
#' @seealso [pairwiseLR()]
#'
#' @return A list with the following entries:
#'
#'   * `undisputed`: A list of undisputed matches and the corresponding LR
#'   values.
#'
#'   * `pmReduced`: Same as `pm`, but with the undisputed victims removed.
#'
#'   * `amReduced`: Same as `am`, but with the data from undisputed victims
#'   inserted for the corresponding missing persons.
#'
#'   * `missingReduced`: Same as `missing`, but without the undisputed
#'   identified missing persons.
#'
#'   * `LRmatrix`, `LRlist`, `pairings`: Output from `pairwiseLR()` applied to
#'   the reduced problem.
#'
#' @examples
#'
#' pm = planecrash$pm
#' am = planecrash$am
#' missing = planecrash$missing
#'
#' findUndisputed(pm, am, missing, threshold = 1e4)
#'
#' # With `relax = TRUE`, one more identification is undisputed
#' findUndisputed(pm, am, missing, threshold = 1e4, relax = TRUE)
#'
#' @export
findUndisputed = function(pm, am, missing, pairings = NULL, ignoreSex = FALSE, threshold = 10000, relax = FALSE, 
                          limit = 0, check = TRUE, verbose = TRUE) {
  
  if(is.singleton(pm))
    pm = list(pm)
  
  # Victim labels
  vics = unlist(labels(pm))
  names(pm) = vics  # ensure pm is named
  
  # Initialise output
  RES = list()
  
  it = 0
  
  # Pairwise LR matrix
  ss = pairwiseLR(pm, am, missing, pairings = pairings, ignoreSex = ignoreSex, check = check, limit = limit)
  B = ss$LRmatrix
  
  # Loop until problem solved - or no more undisputed matches
  while(length(missing) > 0 && length(vics) > 0 && any(B <= threshold)) {
    
    if(verbose)
      message("\nIteration ", it <- it+1, ":")
      
    # Indices of matches exceeding threshold
    highIdx = which(B > threshold, arr.ind = TRUE)
    
    if(!relax) { # undisputed = no others in row/column exceed 1
      goodRows = which(rowSums(B <= 1) == ncol(B) - 1)
      goodCols = which(colSums(B <= 1) == nrow(B) - 1)
      isUndisp = highIdx[, "row"] %in% goodRows & highIdx[, "col"] %in% goodCols
    }
    else { # undisputed = no others in row/column exceed LR/threshold
      isUndisp = sapply(seq_len(nrow(highIdx)), function(k) {  # safer than apply(.., 1)!
        rw = highIdx[k,1]
        cl = highIdx[k,2]
        all(c(B[rw, -cl], B[-rw, cl]) <= B[rw, cl]/threshold)
      })
    }
      
    if(!any(isUndisp)) {
      if(verbose) message("No undisputed matches")
      break
    }
    
    undisp = highIdx[isUndisp, , drop = FALSE]
    
    for(i in seq_len(nrow(undisp))) {
      rw = undisp[i,1]
      cl = undisp[i,2]
      vic = vics[rw] 
      RES[[vic]] = list(match = missing[cl], LR = B[rw,cl])
      
      if(verbose)
        message(sprintf(" %s = %s (LR = %.3g)", vic, missing[cl], B[rw,cl]))
    }
    
    undispVics = vics[undisp[, 1]]
    undispMP = missing[undisp[, 2]]
    
    ### Update the LR matrix
    
    # Move vic data to AM data
    am = transferMarkers(from = pm, to = am, idsFrom = undispVics, idsTo = undispMP, erase = FALSE)
    
    # Remove identified names from vectors
    missing = setdiff(missing, undispMP)
    vics = setdiff(vics, undispVics)
    
    # Remove vic from pm
    pm = pm[vics]
    
    # Update `pairings`, if given
    if(!is.null(pairings))
      pairings = lapply(pairings[vics], function(v) setdiff(v, undispMP))
    
    ss = pairwiseLR(pm, am, missing, pairings = pairings, ignoreSex = ignoreSex, check = FALSE, limit = limit)
    B = ss$LRmatrix
  }
  
  c(list(undisputed = RES, pmReduced = pm, amReduced = am, missingReduced = missing), ss)
}
