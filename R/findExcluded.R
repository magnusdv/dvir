#' Excluded individuals and pairings in a DVI dataset
#'
#' @param dvi A [dviData()] object.
#' @param maxIncomp An integer. A pairing is excluded if the number of
#'   incompatible markers exceeds this.
#' @param pairings (Optional) A list of possible pairings for each victim. By
#'   default, `dvi$pairings` is used, or, if this is NULL,
#'   `generatePairings(dvi, ignoreSex = ignoreSex)`.
#' @param ignoreSex A logical, relevant only if both `pairings` and
#'   `dvi$pairings` are NULL. Default: FALSE.
#' @param verbose A logical, by default TRUE.
#'
#' @return A list with the following entries:
#'
#'   * `exclusionMatrix`: A matrix showing the number of inconsistencies for
#'   each pair (or NA if the pairing was not considered)
#'
#'   * `excluded`: A list of three character vectors:
#'       * `sample`: victim samples excluded against all missing persons
#'       * `missing`: missing persons excluded against all victims
#'       * `fam`: families in which all missing members are excluded against
#'   all victim samples
#'
#'   The above elements are removed in the reduced dataset.
#'
#'   * `dviReduced`: A reduced version of `dvi`, where excluded victims/missing
#'   persons are removed.
#'
#'   * `summary`: A list of data frames `PM` and `AM`, summarising the excluded 
#'   individuals.
#'
#' @seealso [findUndisputed()]
#'
#' @examples
#' findExcluded(icmp)
#'
#' @export
findExcluded = function(dvi, maxIncomp = 2, pairings = NULL, ignoreSex = FALSE, verbose = TRUE) {
  
  if(verbose) {
    cat("Finding exclusions\n")
    cat("Max incompatible markers =", maxIncomp, "\n")
  }
  
  # Ensure proper dviData object
  dvi = consolidateDVI(dvi)
  
  pairings = pairings %||% dvi$pairings %||% generatePairings(dvi, ignoreSex = ignoreSex)
  
  pm = dvi$pm
  am = dvi$am
  missing = dvi$missing
  vics = names(pm)
  

  # Exclusion matrix --------------------------------------------------------

  # Initialise matrix
  mat = matrix(NA_integer_, nrow = length(pm), ncol = length(missing), 
               dimnames = list(vics, missing))
  
  # AM components (for use in output)
  comp = getFamily(dvi, dvi$missing)
  
  # Loop through each pair of victim vs missing
  for(vic in vics) {
    vicdata = pm[[vic]]
    compatMiss = setdiff(pairings[[vic]], "*")
    for(mis in compatMiss) {
      ref = am[[comp[mis]]]
      mat[vic, mis] = length(findExclusions(ref, id = mis, candidate = vicdata))
    }
  }

  # Utilities for analysing rows/cols of `mat`
  minEx = function(v)
    if(all(is.na(v))) Inf else min(v, na.rm = TRUE)
  
  commnt = function(v) 
    ifelse(v == Inf, "no compatible pairings", sprintf("min %s incons", v))
  
  # Victims excluded against all --------------------------------------------
  
  pmMinEx = apply(mat, 1, minEx)
  pmNomatch = pmMinEx > maxIncomp
  excludedVics = vics[pmNomatch]
  
  if(!length(excludedVics)) {
    summaryPM = NULL
    if(verbose)
      cat("\nPM samples excluded against all missing: None\n\n")
  }
  else {
    summaryPM = data.frame(Sample = excludedVics, 
                           Conclusion = "excluded", 
                           Comment = commnt(pmMinEx[pmNomatch]),
                           row.names = NULL)
    if(verbose) {
      cat("\nPM samples excluded against all missing:\n")
      cat(sprintf(" %s (%s)", summaryPM$Sample, summaryPM$Comment), sep = "\n")
      cat("\n")
    }
  }
  
  # Missing excluded against all --------------------------------------------
  
  missMinEx = apply(mat, 2, minEx)
  missNomatch = missMinEx > maxIncomp
  excludedMissing = missing[missNomatch]
  
  if(!length(excludedMissing)) {
    summaryAM = NULL
    if(verbose)
      cat("Missing persons excluded against all PM samples: None\n\n")
  }
  else {
    summaryAM = data.frame(Family = comp[excludedMissing],
                           Missing = excludedMissing,
                           Conclusion = "excluded",
                           Comment = commnt(missMinEx[missNomatch]),
                           row.names = NULL)
    if(verbose) {
      cat("Missing persons excluded against all PM samples:\n")
      cat(sprintf(" %s (%s)", summaryAM$Missing, summaryAM$Comment), sep = "\n")
      cat("\n")
    }
  }
  
    # Reduced DVI problem -----------------------------------------------------
  
  keepVics = vics[!pmNomatch]
  keepMissing = missing[!missNomatch]
  
  # Removed families
  famnames = names(am) %||% 1:length(am)
  excludedFams = setdiff(famnames, comp[keepMissing])
  
  if(length(excludedVics) + length(excludedMissing) > 0) {
    dviRed = subsetDVI(dvi, pm = keepVics, missing = keepMissing, verbose = verbose)
  } else {
    dviRed = dvi
  }
  
  # Updated pairings ----------------------------------------------------
  
  nRemov = sum(mat > maxIncomp, na.rm = TRUE)
  if(nRemov > 0) {
    keepPairs = apply(mat[keepVics, , drop = FALSE], 1, function(rw)
      c("*", missing[!is.na(rw) & rw <= maxIncomp]), simplify = FALSE)
  } else {
    keepPairs = pairings
  }
  if(verbose)
    cat(sprintf("Pairings excluded in total: %d\n", nRemov))
  
  
  dviRed$pairings = keepPairs
    
  # Return list -------------------------------------------------------------

  excluded = list(sample = excludedVics, missing = excludedMissing, fam = excludedFams)
  summary = list(PM = summaryPM, AM = summaryAM)
  
  list(exclusionMatrix = mat, excluded = excluded, dviReduced = dviRed, 
       report = summary, summary = summary)
}



#' Find the number of incompatible markers for each
#'
#' This function computes the number of exclusions, i.e., the number of
#' incompatible markers, for each pairwise comparison. By default, mutation
#' models are ignored. The main work is done by [forrel::findExclusions()].
#' 
#' @param dvi A `dviData` object, typically created with [dviData()].
#' @param removeMut A logical. If TRUE (default), all mutations models are
#'   stripped.
#'
#' @return An integer matrix with `length(pm)` columns and `length(am)` rows.
#'
#' @examples
#'
#' # Plane crash example
#' exclusionMatrix(planecrash)
#'
#' # Inspect a particular pair: M3 vs V6
#' pm = planecrash$pm
#' am = planecrash$am
#' forrel::findExclusions(am, id = "M3", candidate = pm$V6)
#'
#' # Plot one of the incompatible markers
#' plotDVI(planecrash, pm = 6, am = 3, marker ="D7S820")
#'
#' @importFrom forrel findExclusions
#' @export
exclusionMatrix = function(dvi, removeMut = TRUE) {
  
  # Ensure proper dviData object
  dvi = consolidateDVI(dvi)
  
  pm = dvi$pm
  am = dvi$am
  missing = dvi$missing
  
  npm = length(pm)
  nmiss = length(missing)
  
  # Initialise matrix
  mat = matrix(0L, nrow = npm, ncol = nmiss, dimnames = list(names(pm), missing))
  
  # AM components (for use in output)
  comp = getFamily(dvi, dvi$missing)
  
  # Loop through each pair of victim vs missing
  for(i in 1:npm) for(j in 1:nmiss) {
    vic = pm[[i]]
    ref = am[[comp[j]]]
    mat[i, j] = length(findExclusions(ref, id = missing[j], candidate = vic))
  }
  
  mat
}

