#' Undisputed identifications in a DVI problem
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
#' @param dvi A `dviData` object, typically created with [dviData()].
#' @param pairings A list of possible pairings for each victim. If NULL, all
#'   sex-consistent pairings are used.
#' @param ignoreSex A logical.
#' @param threshold A non-negative number. If no pairwise LR exceed this, the
#'   iteration stops.
#' @param relax A logical affecting the definition of being undisputed (see
#'   Details). Default: FALSE.
#' @param limit A positive number. Only pairwise LR values above this are
#'   considered.
#' @param nkeep An integer, or NULL. If given, only the `nkeep` most likely
#'   pairings are kept for each victim.
#' @param check A logical indicating if the input data should be checked for
#'   consistency. Default: TRUE.
#' @param numCores An integer; the number of cores used in parallelisation.
#'   Default: 1.
#' @param verbose A logical. Default: TRUE.
#'
#' @seealso [pairwiseLR()], [findExcluded()]
#'
#' @return A list with the following entries:
#'
#'   * `undisputed`: A data frame containing the undisputed matches, including LR.
#'
#'   * `dviReduced`: A reduced version of `dvi`, where undisputed
#'   victims/missing persons are removed, and data from undisputed victims
#'   inserted in `am`.
#'
#'   * `LRmatrix`, `LRlist`, `pairings`: Output from `pairwiseLR()` applied to
#'   the reduced problem.
#'
#' @examples
#'
#' \donttest{
#' findUndisputed(planecrash, threshold = 1e4)
#'
#' # With `relax = TRUE`, one more identification is undisputed
#' findUndisputed(planecrash, threshold = 1e4, relax = TRUE)
#' }
#'
#' @export
findUndisputed = function(dvi, pairings = NULL, ignoreSex = FALSE, threshold = 10000, 
                          relax = FALSE, limit = 0, nkeep = NULL, check = TRUE, 
                          numCores = 1, verbose = TRUE) {
  
  if(verbose) {
    cat("\nFinding undisputed matches\n")
    cat("Pairwise LR threshold =", threshold, "\n")
  }

  # Ensure proper dviData object
  dvi = consolidateDVI(dvi)
  
  # AM components (for use in output)
  comp = getFamily(dvi, dvi$missing)
  
  # Initialise pairings if not given
  dvi$pairings = pairings %||% dvi$pairings %||% generatePairings(dvi, ignoreSex = ignoreSex)
  
  # Initialise output
  RES = list()
  
  # Loop until problem solved - or no more undisputed matches
  stp = 0
  while(TRUE) {
    
    if(verbose)
      cat(sprintf("\nStep %d:\n", stp <- stp+1))
    
    vics = names(dvi$pm)
    missing = dvi$missing
    
    # Pairwise LR matrix
    ss = pairwiseLR(dvi, check = check, limit = limit, nkeep = nkeep, 
                    numCores = numCores, verbose = verbose)
    B = ss$LRmatrix
    
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
     
    Nundisp = if(!length(isUndisp)) 0 else sum(isUndisp)
    
    if(!Nundisp) {
      if(verbose)
        cat(sprintf("No%s undisputed matches\n", if(stp > 1) " further" else ""))
      break
    }
    
    if(verbose)
      cat(sprintf("%d undisputed match%s\n", Nundisp, if(Nundisp == 1) "" else "es"))
    
    undisp = highIdx[isUndisp, , drop = FALSE]
    
    for(i in seq_len(Nundisp)) {
      rw = undisp[i,1]
      cl = undisp[i,2]
      vic = vics[rw] 
      RES[[vic]] = list(Step = stp, Missing = missing[cl], LR = B[rw,cl])
      
      if(verbose)
        cat(sprintf(" %s = %s (LR = %.3g)\n", vic, missing[cl], B[rw,cl]))
    }
    
    undispVics = vics[undisp[, 1]]
    undispMP = missing[undisp[, 2]]
    
    # Data from identified samples (keep for updating dviRed below)
    undispData = dvi$pm[undispVics]
    
    # Reduced DVI dataset
    newvics = setdiff(vics, undispVics)
    newmissing = setdiff(missing, undispMP)
    dvi = subsetDVI(dvi, pm = newvics, missing = newmissing, verbose = verbose)
    
    # Move vic data to AM data
    names(undispVics) = undispMP
    relevantMP = intersect(undispMP, unlist(labels(dvi$am)))
    if(length(relevantMP))
      dvi$am = transferMarkers(from = undispData, to = dvi$am, 
                               idsFrom = undispVics[relevantMP], 
                               idsTo = relevantMP, erase = FALSE)
    
    # Break?
    if(!length(newmissing) || !length(newvics))
      break
  }
  
  undispDF = do.call(rbind.data.frame, RES)
  if(nrow(undispDF)) {
    undispDF = cbind(undispDF, Sample = names(RES), Family = comp[undispDF$Missing])
    undispDF = undispDF[c("Sample", "Missing", "Family", "LR", "Step")]
    rownames(undispDF) = NULL
  }
  
  # Disputed: TODO!
  
  c(list(undisputed = undispDF, dviReduced = dvi), ss)
}
