#' Sequential DVI search
#'
#' Two sequential approaches based on the marginal LR's are implemented. See
#' Details for explanations.
#'
#' Both methods proceed iteratively, always choosing the match `V_i = M_j` with
#' highest marginal LR. If the maximum is not unique, the first match is chosen.
#' The two methods differ only in how the LR matrix is updated in each step.
#'
#' In `sequential1()`, the same marginal matrix is used in all steps, with the
#' only modification that previously matching individuals are ignored. (After
#' each identification, the corresponding row and column are set to 0.)
#'
#' In `sequential2()` the marginal LR matrix is re-computed in each step, after
#' transferring the identified victim into the AM data.
#'
#' @param pm PM data: List of singletons.
#' @param am AM data: A ped object or list of such.
#' @param MPs Character vector with names of the missing persons.
#' @param threshold A non-negative number. If no marginal LR values exceed this,
#'   the iteration stops.
#' @param check A logical, indicating if the input data should be checked for consistency.
#' @param verbose A logical.
#'
#' @return A character vector of the same length as `pm`, indicating the
#'   solution. For example, the result vector `c("MP2", "*", "MP3)` corresponds
#'   to the assignment `V1 = MP2, V3 = MP3`.
#'
#' @examples
#'
#' pm = example1$pm
#' am = example1$am
#' MPs = example1$MPs
#'
#' sequential1(pm, am, MPs)
#'
#' sequential2(pm, am, MPs)
#'
#'
#' @export
sequential1 = function(pm, am, MPs, threshold = 1, check = TRUE, verbose = FALSE) {
  # Victim labels
  vics = unlist(labels(pm))
  
  # Marginal matrix
  marg = marginal(pm, am, MPs, check = check)$LR.table
  
  # Initialise solution vector with no moves
  RES = rep("*", length(pm))
  names(RES) = vics
  
  i = 0
  
  # Loop until all LRs are below threshold or all victims are identified
  while(any(marg > threshold)) {
    
    # Find cell with highest LR (if not unique, take the first!)
    mx = which(marg == max(marg), arr.ind = TRUE)[1, ]
    vic = vics[mx[1]]
    mp = MPs[mx[2]]
    
    RES[vic] = mp
    
    if(verbose) {
      cat("Iteration", i <- i+1, "\n")
      print(marg)
      message(sprintf("--> %s = %s", vic, mp))
    }
    
    # Mute corresponding row & column
    marg[vic, ] = marg[, mp] = 0
  }
  
  RES
}


#' @rdname sequential1
#' @export
sequential2 = function(pm, am, MPs, threshold = 1, check = TRUE, verbose = FALSE) {
  # Victim labels
  vics = unlist(labels(pm))
  
  # Initialise solution vector with no moves
  RES = rep("*", length(pm))
  names(RES) = vics
  
  i = 0
  
  # Loop until all LRs are below threshold or all victims are identified
  while(length(MPs) > 0) {
    # Marginal matrix
    marg = marginal(pm, am, MPs, check = check)$LR.table
    
    # If no matches: stop
    if(all(marg < threshold))
      break
    
    # Find cell with highest LR (if not unique, take the first!)
    mx = which(marg == max(marg), arr.ind = TRUE)[1, ]
    vic = vics[mx[1]]
    mp = MPs[mx[2]]
    
    RES[vic] = mp
    
    if(verbose) {
      cat("Iteration", i<-i+1, "\n")
      print(marg)
      message(sprintf("--> %s = %s", vic, mp))
    }
    
    ### Update the marginal LR matrix
    
    # Move vic data to AM data
    am = transferMarkers(from = pm, to = am, idsFrom = vic, idsTo = mp)

    # Remove identified names from vectors
    MPs = setdiff(MPs, mp)
    vics = setdiff(vics, vic)
    
    # Remove vic from pm
    pm = pm[vics]
  }
  
  RES
}