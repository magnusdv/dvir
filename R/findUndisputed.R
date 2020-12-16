#' Undisputed identifications in DVI problems
#'
#' This function uses the single-search LR values to find "undisputed" matches
#' between victims and missing individuals. An identification \eqn{V_i = M_j} is
#' called undisputed if the corresponding likelihood ratio \eqn{LR_{i,j}}
#' exceeds the given threshold, while all other values involving \eqn{v_i} or
#' \eqn{M_j} are below 1.
#'
#' @param pm PM data: List of singletons.
#' @param am AM data: A `ped` object or list of such.
#' @param missing Character vector with names of the missing persons.
#' @param threshold A non-negative number. If no single-search LR values exceed
#'   this, the iteration stops.
#' @param limit A positive number. Only single-search LR values above this are
#'   considered.
#' @param check A logical, indicating if the input data should be checked for
#'   consistency.
#' @param verbose A logical.
#'
#' @return A list with the following entries:
#'
#'   * `undisputed`: A list of undisputed matches and the corresponding LR
#'   values.
#'
#'   * `pmReduced`: Same as `pm`, but with the undisputed victims removed.
#'
#'   * `amReduced`: Same as `am`, but with the data from undisputed victims
#'   inserted.
#'
#'   * `missingReduced`: Same as `missing`, but without the undisputed
#'   identified missing persons.
#'
#'   * `moves`, `LR.table`, `LRmatrix`: Output from `singleSearch()` applied to
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
#' @export
findUndisputed = function(pm, am, missing, threshold = 10000, limit = 0, check = TRUE, verbose = FALSE) {
  
  if(is.singleton(pm))
    pm = list(pm)
  
  # Victim labels
  vics = unlist(labels(pm))
  names(pm) = vics  # ensure pm is named
  
  # Initialise output
  RES = list()
  
  it = 0
  
  # single-search matrix
  ss = singleSearch(pm, am, missing, check = check)
  marg = ss$LR.table
  
  # Loop until problem solved - or no more undisputed matches
  while(length(missing) > 0 && length(vics) > 0 && any(marg <= threshold)) {
    
    if(verbose)
      message("\nIteration ", it <- it+1, ":")
      
    # Indices of matches exceeding threshold
    highIdx = which(highLR <- marg > threshold, arr.ind = TRUE)
    
    # Find "undisputed" matches, i.e., no others in row/column exceed 1
    goodRows = which(rowSums(marg <= 1) == ncol(marg) - 1)
    goodCols = which(colSums(marg <= 1) == nrow(marg) - 1)
    isUndisp = highIdx[, "row"] %in% goodRows & highIdx[, "col"] %in% goodCols
    
    if(!any(isUndisp)) {
      if(verbose) message("No undisputed matches")
      break
    }
    
    undisp = highIdx[isUndisp, , drop = FALSE]
    
    for(i in seq_len(nrow(undisp))) {
      rw = undisp[i,1]
      cl = undisp[i,2]
      vic = vics[rw] 
      RES[[vic]] = list(match = missing[cl], LR = marg[rw,cl])
      
      if(verbose)
        message(sprintf(" %s = %s (LR = %.3g)", vic, missing[cl], marg[rw,cl]))
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
    
    ss = singleSearch(pm, am, missing, check = FALSE)
    marg = ss$LR.table
  }
  
  c(list(undisputed = RES, pmReduced = pm, amReduced = am, missingReduced = missing), ss)
}
