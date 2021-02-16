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
    method = sprintf("Method: Sequential search %s LR updates\n", ifelse(updateLR, "with", "without"))
    summariseDVI(pm, am, missing, printMax = 10, method = method)
  }
  
  # Initialise solution vector with no moves
  sol = rep("*", length(pm))
  
  # Ensure pm and sol is properly named
  names(sol) = names(pm) = unlist(labels(pm)) 
  
  # LR matrix
  B = pairwiseLR(pm, am, missing, check = check)$LRmatrix
  
  # Environment for keeping parameters and storing solutions
  env = list2env(list(RES = list(), updateLR = updateLR, threshold = threshold, 
                      verbose = verbose, debug = debug))
  
  # Start recursion
  addPairing(pm, am, B, sol, env)
  
  as.data.frame(do.call(rbind, unique(env$RES)))
}

# Recursive function, adding one new pairing to the solution vector `sol`
# If final: Store in the RES list of the environment `env`
addPairing = function(pm, am, B, sol, env) {
  
  step = length(sol) - nrow(B)
  
  # If all below threshold: store current solution and stop
  if(all(B < env$threshold - sqrt(.Machine$double.eps))) {
    env$RES = append(env$RES, list(sol))
    if(env$verbose) 
      message(sprintf("%sStep %d: Stop (max LR = %.2f)", strrep(" ", step), step + 1, max(B)))
    if(env$debug) {print(B); cat("\n")}
    return()
  }
  
  vics = rownames(B)
  missing = colnames(B)
  
  # Indices of maximal elements
  allmax = which(B == max(B), arr.ind = TRUE)
  
  # For each max, set as solution and recurse
  # If multiple max'es, this creates new branches
  for(i in 1:nrow(allmax)) {
    mx = allmax[i, ]
    vic = vics[mx[1]]
    mp = missing[mx[2]]
    sol[vic] = mp
    
    if(env$verbose) {
      message(sprintf("%sStep %d: %s = %s (LR = %.2g)", 
                      strrep(" ", step), step + 1, vic, mp, max(B)))
      if(env$debug) {print(B); cat("\n")}
    }
    
    # If all victims identified: store solution and return
    if(min(dim(B)) == 1) {
      env$RES = append(env$RES, list(sol))
      return()
    }
    
    ### Update DVI: Transfer vic data from old PM to new AM
    newPm = pm[-mx[1]]
    newAm = transferMarkers(from = pm, to = am, idsFrom = vic, idsTo = mp, erase = FALSE)
    newMissing = setdiff(missing, mp)
    
    if(env$updateLR)
      newB = pairwiseLR(newPm, newAm, newMissing, check = FALSE)$LRmatrix
    else
      newB = B[-mx[1], -mx[2], drop = FALSE]
    
    # Recurse
    addPairing(newPm, newAm, newB, sol, env)
  }
}