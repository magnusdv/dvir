#' Sequential DVI search
#'
#' @param pm PM data: List of singletons.
#' @param am AM data: A ped object or list of such.
#' @param missing Character vector with names of the missing persons.
#' @param updateLR A logical. If TRUE, the LR matrix is updated in each
#'   iteration.
#' @param threshold A non-negative number. If no pairwise LR values exceed
#'   this, the iteration stops.
#' @param check A logical, indicating if the input data should be checked for
#'   consistency.
#' @param verbose A logical.
#' @param debug A logical. If TRUE, the LR matrix is printed 
#' 
#'
#' @return A solution to the DVI problem in the form of an assignment vector.
#'
#' @examples
#' pm = example1$pm
#' am = example1$am
#' missing = example1$missing
#'
#' sequentialDVI(pm, am, missing, updateLR = FALSE)
#' sequentialDVI(pm, am, missing, updateLR = TRUE)
#'
#' @export
sequentialDVI = function(pm, am, missing, updateLR = TRUE, 
                         threshold = 1, check = TRUE, verbose = TRUE, debug = FALSE) {
  
  if(is.singleton(pm))
    pm = list(pm)
  
  if(verbose) {
    method = sprintf("Sequential search %s LR updates", ifelse(updateLR, "with", "without"))
    summariseDVI(pm, am, missing, printMax = 10, method = method)
  }
  
  # Victim labels
  names(pm) = vics = unlist(labels(pm)) 
  
  # Initialise solution vector with no moves
  RES = rep("*", length(pm))
  names(RES) = vics
  
  # LR matrix
  marg = pairwiseLR(pm, am, missing, check = check)$LRmatrix
  
  i = 0
  
  # Loop until all LRs are below threshold or all victims are identified
  while(length(missing) > 0) {
    
    # If no matches: stop
    if(all(marg < threshold))
      break
    
    # Find cell with highest LR (if not unique, pick random)
    mx = which(marg == max(marg), arr.ind = TRUE)
    if(nrow(mx) > 1)
      mx = mx[sample(nrow(mx), size = 1), ]
    
    vic = vics[mx[1]]
    mp = missing[mx[2]]
    
    RES[vic] = mp
    
    if(verbose) {
      message("\nIteration ", i<-i+1)
      if(debug) print(marg)
      message(sprintf("---> %s = %s (LR = %.2g)", vic, mp, max(marg)))
    }
    
    if(i == length(RES)) 
      break 
    
    ### Update the LR matrix
    
    if(updateLR) {
      
      # Move vic data to AM data
      am = transferMarkers(from = pm, to = am, idsFrom = vic, idsTo = mp, erase = FALSE)
      
      # Remove identified names from vectors
      missing = setdiff(missing, mp)
      vics = setdiff(vics, vic)
      
      # Remove vic from pm
      pm = pm[vics]
      
      # Re-compute the LR matrix
      marg = pairwiseLR(pm, am, missing, check = check)$LRmatrix
    }
    else {
      # Mute corresponding row & column
      marg[vic, ] = marg[, mp] = 0
    }
  }
  
  RES
}
