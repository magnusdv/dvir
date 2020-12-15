#' Sequential DVI search
#'
#' @param pm PM data: List of singletons.
#' @param am AM data: A ped object or list of such.
#' @param missing Character vector with names of the missing persons.
#' @param updateLR A logical. If TRUE, the LR matrix is updated in each
#'   iteration.
#' @param threshold A non-negative number. If no single-search LR values exceed
#'   this, the iteration stops.
#' @param check A logical, indicating if the input data should be checked for
#'   consistency.
#' @param verbose A logical.
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
                         threshold = 1, check = TRUE, verbose = FALSE) {
  
  if(is.singleton(pm))
    pm = list(pm)
  
  # Victim labels
  vics = unlist(labels(pm))
  
  # Initialise solution vector with no moves
  RES = rep("*", length(pm))
  names(RES) = vics
  
  # LR matrix
  marg = singleSearch(pm, am, missing, check = check)$LR.table
  
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
      cat("Iteration", i<-i+1, "\n")
      print(marg)
      message(sprintf("--> %s = %s", vic, mp))
    }
    
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
      marg = singleSearch(pm, am, missing, check = check)$LR.table
    }
    else {
      # Mute corresponding row & column
      marg[vic, ] = marg[, mp] = 0
    }
  }
  
  RES
}
