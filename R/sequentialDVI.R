#' Sequential DVI search
#'
#' Performs a sequential matching procedure based on the pairwise LR matrix. In
#' each step the pairing corresponding to the highest LR is selected and
#' included as a match if the LR exceeds the given threshold. By default,
#' (`updateLR = TRUE`) the pairwise LRs are recomputed in each step after
#' including the data from the identified sample.
#'
#' If, at any point, the highest LR is obtained by more than one pairing, the
#' process branches off and produces multiple solutions. (See Value.)
#'
#' @param dvi A `dviData` object, typically created with [dviData()].
#' @param updateLR A logical. If TRUE, the LR matrix is updated in each
#'   iteration.
#' @param threshold A non-negative number. If no pairwise LR values exceed this,
#'   the iteration stops.
#' @param check A logical, indicating if the input data should be checked for
#'   consistency.
#' @param verbose A logical.
#' @param debug A logical. If TRUE, the LR matrix is printed.
#'
#' @return A list with two elements:
#'  * `matches`: A single assignment vector, or (if multiple branches)
#'   a data frame where each row is an assignment vector.
#'  * `details`: A data frame (or a list of data frames, if multiple branches)
#'   including the LR of each identification in the order they were made.
#'
#' @examples
#' # Without LR updates
#' sequentialDVI(example1, updateLR = FALSE)
#'
#' # With LR updates (default). Note two branches!
#' r = sequentialDVI(example1)
#'
#' # Plot the two solutions
#' plotSolution(example1, r$matches, k = 1)
#' plotSolution(example1, r$matches, k = 2)
#' 
#' # Add `debug = T` to see the LR matrix in each step
#' sequentialDVI(example1, debug = TRUE)
#' 
#' # The output of can be fed into `jointDVI()`:
#' jointDVI(example1, assignments = r$matches)
#' 
#' @importFrom stats setNames
#' @export
sequentialDVI = function(dvi, updateLR = TRUE, threshold = 1, check = TRUE, 
                         verbose = TRUE, debug = FALSE) {
  
  # Ensure proper dviData object
  dvi = consolidateDVI(dvi)
  vics = names(dvi$pm)
  nVics = length(vics)
  
  # AM components (for use in output)
  comp = getFamily(dvi, dvi$missing)
  
  if(verbose) {
    print.dviData(dvi)
    cat(sprintf("\nSequential search %s LR updates\n", ifelse(updateLR, "with", "without")))
  }
  
  # Initialise 'null' solution
  matches = list() 
  
  # LR matrix
  B = pairwiseLR(dvi, check = check)$LRmatrix
  
  # Environment for keeping parameters and storing solutions
  env = list2env(list(RES = list(), updateLR = updateLR, threshold = threshold, 
                      nVics = nVics, verbose = verbose, debug = debug))
  
  # Start recursion
  addPairing(dvi, B, matches, env)
  
  # Collect result table(s) with LR
  restabs = lapply(env$RES, function(lst) {
    df = do.call(rbind.data.frame, lst)
    if(nrow(df)) {
      df = cbind(df, Sample = names(lst), Family = comp[df$Missing])
      df = df[order(df$Step), c("Sample", "Missing", "Family", "LR", "Step")]
      rownames(df) = NULL
    }
    df
  })
  
  # Extract assignments
  A = setNames(rep("*", nVics), vics)
  mList = lapply(restabs, function(tab) {a = A; a[tab$Sample] = tab$Missing; a})
  matches = as.data.frame(do.call(rbind, unique(mList)))
  
  res = list(matches = matches, details = restabs)
  
  # Simplify output if no branching
  if(length(restabs) == 1)
    res = list(matches = matches[1,], details = restabs[[1]])
  
  res
}

# Recursive function, adding one new identification to the `matches` list
# If final: Store in the RES list of the environment `env`
addPairing = function(dvi, B, matches, env) {
  step = env$nVics - nrow(B)
  Bmax = max(B)
  
  # If all below threshold: store current solution and stop
  if(Bmax < env$threshold - sqrt(.Machine$double.eps)) {
    env$RES = append(env$RES, list(matches))
    if(env$verbose) 
      cat(sprintf("%sStep %d: Stop (max LR = %.2f)\n", strrep(" ", step), step + 1, max(B)))
    if(env$debug) {print(B); cat("\n")}
    return()
  }
  
  vics = rownames(B)
  missing = colnames(B)
  
  # Indices of maximal elements
  allmax = which(B == Bmax, arr.ind = TRUE)
  
  # For each max, set as solution and recurse
  # If multiple max'es, this creates new branches
  for(i in 1:nrow(allmax)) {
    mx = allmax[i, ]
    vic = vics[mx[1]]
    mp = missing[mx[2]]
    matches[[vic]] = list(Missing = mp, LR = Bmax, Step = step+1)
    if(env$verbose) {
      cat(sprintf("%sStep %d: %s = %s (LR = %.2g)\n", 
                      strrep(" ", step), step + 1, vic, mp, Bmax))
      if(env$debug) {print(B); cat("\n")}
    }
    
    # If all victims identified: store solution and return
    if(min(dim(B)) == 1) {
      env$RES = append(env$RES, list(matches))
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
    addPairing(newDvi, newB, matches, env)
  }
}
