#' Sequential DVI search
#'
#' @param dvi A `dviData` object, typically created with [dviData()].
#' @param updateLR A logical. If TRUE, the LR matrix is updated in each
#'   iteration.
#' @param threshold A non-negative number. If no pairwise LR values exceed this,
#'   the iteration stops.
#' @param check A logical, indicating if the input data should be checked for
#'   consistency.
#' @param verbose A logical.
#' @param debug A logical. If TRUE, the LR matrix is printed
#'
#'
#' @return A solution to the DVI problem in the form of an assignment vector.
#'
#' @examples
#' sequentialDVI(example1, updateLR = FALSE)
#' sequentialDVI(example1, updateLR = TRUE)
#'
#' # The output of can be fed into `jointDVI()`:
#' res = sequentialDVI(example1, updateLR = TRUE)
#' jointDVI(example1, assignments = res)
#' 
#' @export
sequentialDVI = function(dvi, updateLR = TRUE, threshold = 1, 
                         check = TRUE, verbose = TRUE, debug = FALSE) {
  
  if(!inherits(dvi, "dviData"))
    stop2("First argument must be `dviData` object. (As of dvir version 3.0.0)")
  
  if(is.singleton(dvi$pm))
    dvi$pm = list(dvi$pm)
  
  if(verbose) {
    method = sprintf("Method: Sequential search %s LR updates\n", ifelse(updateLR, "with", "without"))
    summariseDVI(dvi, printMax = 10, method = method)
  }
  
  # Initialise 'null' solution
  sol = rep("*", length(dvi$pm))
  
  # Ensure pm and sol is properly named
  names(sol) = names(dvi$pm) = unlist(labels(dvi$pm)) 
  
  # LR matrix
  B = pairwiseLR(dvi, check = check)$LRmatrix
  
  # Environment for keeping parameters and storing solutions
  env = list2env(list(RES = list(), updateLR = updateLR, threshold = threshold, 
                      verbose = verbose, debug = debug))
  
  # Start recursion
  addPairing(dvi, B, sol, env)
  
  # Final output
  as.data.frame(do.call(rbind, unique(env$RES)))
}

# Recursive function, adding one new pairing to the solution vector `sol`
# If final: Store in the RES list of the environment `env`
addPairing = function(dvi, B, sol, env) {
  
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
    pm = dvi$pm
    am = dvi$am
    newPm = pm[-mx[1]]
    newAm = transferMarkers(from = pm, to = am, idsFrom = vic, idsTo = mp, erase = FALSE)
    newMissing = setdiff(missing, mp)
    newDvi = dviData(newPm, newAm, newMissing)
    
    if(env$updateLR)
      newB = pairwiseLR(newDvi, check = FALSE)$LRmatrix
    else
      newB = B[-mx[1], -mx[2], drop = FALSE]
    
    # Recurse
    addPairing(newDvi, newB, sol, env)
  }
}
