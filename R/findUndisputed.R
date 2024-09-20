#' Undisputed identifications in a DVI problem
#'
#' This function uses the pairwise LR matrix to find *undisputed* matches
#' between victims and missing individuals. An identification \eqn{V_i = M_j} is
#' called undisputed, relative to a threshold T, if the corresponding likelihood
#' ratio \eqn{LR_{i,j} \geq T} AND \eqn{LR_{i,j}} is at least T times greater
#' than all other pairwise LRs involving \eqn{V_i} or \eqn{M_j}.
#'
#' If the parameter `strict` is set to TRUE, the last criterion is replaced with
#' the stronger requirement that all other pairwise LRs involving \eqn{V_i} or
#' \eqn{M_j} must be at most 1.
#'
#' @param dvi A `dviData` object, typically created with [dviData()].
#' @param pairings A list of possible pairings for each victim. If NULL, all
#'   sex-consistent pairings are used.
#' @param ignoreSex A logical.
#' @param threshold A non-negative number. If no pairwise LR exceed this, the
#'   iteration stops.
#' @param strict A logical affecting the definition of being undisputed (see
#'   Details). Default: FALSE.
#' @param relax Deprecated; use `strict = FALSE` instead.
#' @param limit A positive number. Only pairwise LR values above this are
#'   considered.
#' @param nkeep An integer, or NULL. If given, only the `nkeep` most likely
#'   pairings are kept for each victim.
#' @param numCores An integer; the number of cores used in parallelisation.
#'   Default: 1.
#' @param verbose A logical. Default: TRUE.
#'
#' @seealso [pairwiseLR()], [findExcluded()]
#'
#' @return A list with the following entries:
#'
#'  * `dviReduced`: A reduced version of `dvi`, where undisputed
#'   victims/missing persons are removed, and data from undisputed victims
#'   inserted into the reference data.
#'
#'  * `summary`: A data frame summarising the undisputed matches.
#'
#'  * `LRmatrix`: Output from `pairwiseLR()` applied to
#'   the reduced problem.
#'
#' @examples
#'
#' \donttest{
#' u1 = findUndisputed(planecrash, verbose = FALSE)
#' u1$summary 
#' 
#' # With `strict = TRUE`, the match M3 = V2 goes away
#' u2 = findUndisputed(planecrash, strict = TRUE, verbose = FALSE)
#' u2$summary
#' 
#' # Reason: M3 has LR > 1 also against V7
#' u2$LRmatrix[, "M3"] |> round(2)
#' }
#'
#' @export
findUndisputed = function(dvi, pairings = NULL, ignoreSex = FALSE, 
                          threshold = 10000, strict = FALSE, relax = !strict, 
                          limit = 0, nkeep = NULL, numCores = 1, verbose = TRUE) {
  
  if(!missing(relax)) {
    cat("Warning: `relax` is deprecated; replaced by (its negation) `strict`")
    strict = !relax
  }
  
  if(verbose) {
    cat("Finding undisputed matches\n")
    cat("Pairwise LR threshold =", threshold, "\n")
  }

  if(!isTRUE(length(threshold) == 1 && threshold >= 0))
    stop2("`treshold` must be a nonnegative number: ", threshold %||% "NULL")
  
  # Ensure proper dviData object
  dvi = consolidateDVI(dvi)
  
  if(!length(dvi$pm) || !length(dvi$missing))
    return(list(dviReduced = dvi, summary = NULL))
  
  # AM components (for use in output)
  comp = getFamily(dvi, dvi$missing)
  
  # Initialise pairings if not given
  dvi$pairings = pairings %||% dvi$pairings %||% generatePairings(dvi, ignoreSex = ignoreSex)
  
  # Initialise output
  RES = list()
  
  # Loop until problem solved - or no more undisputed matches
  stp = 0
  while(TRUE) {
    
    stp = stp + 1
    
    if(verbose)
      cat(sprintf("\nStep %d:\n", stp))
    
    vics = names(dvi$pm)
    missing = dvi$missing
    
    # Pairwise LR matrix
    ss = pairwiseLR(dvi, limit = limit, nkeep = nkeep, numCores = numCores, 
                    check = FALSE, verbose = verbose)
    B = ss$LRmatrix
    
    # Indices of undisputed matches
    undisp = undisputedEntries(B, threshold, strict = strict)
    nu = nrow(undisp)
    
    if(verbose) {
      if(nu == 0) 
        cat(sprintf("No%s undisputed matches\n", if(stp > 1) " further" else ""))
      else 
        cat(sprintf("%d undisputed match%s\n", nu, if(nu == 1) "" else "es"))
    }
    
    if(nu == 0)
      break
    
    for(i in seq_len(nu)) {
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
    
    # Victims removed with no remaining pairings
    for(v in setdiff(newvics, names(dvi$pm)))
      RES[[v]] = list(Step = stp, Missing = NA, LR = NA)
    
      
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
  
  summary = do.call(rbind.data.frame, RES)
  isExcl = is.na(summary$Missing)
  step = paste("Step", summary$Step)
  
  if(nrow(summary)) {
    summary$Sample = names(RES)
    summary$Family = comp[summary$Missing]
    summary$Conclusion = ifelse(isExcl, "Excluded", "Undisputed")
    summary$Comment = ifelse(isExcl, "No compatible pairings", step)
    summary = summary[c("Family", "Missing", "Sample", "LR", "Conclusion", "Comment")]
    rownames(summary) = NULL
  }
  
  # Update pairings using output from the last pairwiseLR
  dvi$pairings = ss$pairings
  
  list(dviReduced = dvi, summary = summary, LRmatrix = B)
}
