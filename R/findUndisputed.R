#' Undisputed identifications in DVI problems
#'
#' This function uses the single-search LR values to find "undisputed" matches
#' between victims and missing individuals. An identification `V_i = M_j` is
#' called undisputed if the corresponding likelihood ratio `LR_{i,j}` exceeds
#' the given threshold, while all other values involving `v_i` or `M_j` are
#' below 1.
#' @param pm PM data: List of singletons.
#' @param am AM data: A `ped` object or list of such.
#' @param missing Character vector with names of the missing persons.
#' @param threshold A non-negative number. If no marginal LR values exceed this,
#'   the iteration stops.
#' @param check A logical, indicating if the input data should be checked for
#'   consistency.
#' @param verbose A logical.
#'
#' @return A named character vector.
#'
#' @examples
#' 
#' pm = planecrash$pm
#' am = planecrash$am
#' missing = planecrash$MPs
#'
#' findUndisputed(pm, am, missing, threshold = 1e4)
#'
#' @export
findUndisputed = function(pm, am, missing, threshold = 10000, check = TRUE, verbose = FALSE) {
  
  if(is.singleton(pm))
    pm = list(pm)
  
  # Victim labels
  vics = unlist(labels(pm))
  
  # Store these for later
  origPM = pm
  origAM = am
  origMissing = missing
  origVics = vics
  
  # Initialise solution vector with no moves
  RES = rep("*", length(pm))
  names(RES) = vics
  
  # Ensure pm is named
  names(pm) = vics
  
  ### Step 1: Sequential
  if(verbose)
    message("\nSequential identification of undisputed matches")
  it = 0
  
  # Loop until problem solved - or no more undisputed matches
  while(length(missing) > 0 && length(vics) > 0) {
    
    # Marginal matrix
    marg = marginal(pm, am, missing, check = check)$LR.table
    if(all(marg <= threshold))
      break
    
    # Indices of matches exceeding threshold
    highIdx = which(highLR <- marg > threshold, arr.ind = TRUE)
    
    # Find "undisputed" matches, i.e., no others in row/column exceed 1
    goodRows = which(rowSums(marg <= 1) == ncol(marg) - 1)
    goodCols = which(colSums(marg <= 1) == nrow(marg) - 1)
    isUndisp = highIdx[, "row"] %in% goodRows & highIdx[, "col"] %in% goodCols
    
    if(!any(isUndisp)) 
      break
    
    undisp = highIdx[isUndisp, , drop = FALSE]
    
    if(verbose) {
      message("\nIteration ", it <- it+1)
      print(marg)
      message("\nUndisputed matches in this iteration:")
      for(i in seq_len(nrow(undisp))) {
        rw = undisp[i,1]; cl = undisp[i,2]
        message(sprintf(" %s = %s (LR = %.3g)", vics[rw], missing[cl], marg[rw,cl]))
      }
    }
    
    undispVics = vics[undisp[, 1]]
    undispMP = missing[undisp[, 2]]
    
    RES[undispVics] = undispMP
    
    ### Update the marginal LR matrix
    
    # Move vic data to AM data
    am = transferMarkers(from = pm, to = am, idsFrom = undispVics, idsTo = undispMP, erase = FALSE)
    
    # Remove identified names from vectors
    missing = setdiff(missing, undispMP)
    vics = setdiff(vics, undispVics)
    
    # Remove vic from pm
    pm = pm[vics]
  }
  
  RES[RES != "*"]
}
