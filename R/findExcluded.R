#' Excluded individuals and pairings in a DVI dataset
#'
#' Analysing exclusions is often an efficient way to reduce large DVI datasets.
#' A pairing V = M is *excluded* if it implies (too many) genetic
#' inconsistencies. The function `findExcluded()` identifies and removes (i)
#' victim samples with too many inconsistencies against all missing persons,
#' (ii) missing persons with too many inconsistencies against all victim
#' samples, and (iii) inconsistent pairings among the remaining.
#'
#' The main calculation in `findExcluded()` is done by `exclusionMatrix()`,
#' which records number of incompatible markers of each pairwise comparison.
#'
#' @param dvi A [dviData()] object.
#' @param maxIncomp An integer. A pairing is excluded if the number of
#'   incompatible markers exceeds this.
#' @param pairings A list of possible pairings for each victim. By default,
#'   `dvi$pairings` is used, or, if this is NULL, `generatePairings(dvi,
#'   ignoreSex)`.
#' @param ignoreSex A logical, by default: FALSE.
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
#'   * `dviReduced`: A reduced version of `dvi`, where the excluded elements
#'   are removed, and the pairings are updated.
#'
#'   * `summary`: A list of data frames `PM` and `AM`, summarising the excluded
#'   individuals.
#'
#' @seealso [findUndisputed()]. See also [forrel::findExclusions()] for analysis
#'   of a specific pairwise comparison.
#'
#' @examples
#' \donttest{
#' e = findExcluded(icmp)
#' e$summary
#' e$exclusionMatrix
#'
#' # The exclusion matrix can also be computed directly:
#' exclusionMatrix(icmp)
#'
#' # Inspect a particular pair: M4 vs V4
#' forrel::findExclusions(icmp$am, id = "M4", candidate = icmp$pm$V4)
#'
#' # Plot one of the incompatible markers
#' plotDVI(icmp, pm = 4, marker ="D7S820")
#' }
#' @export
findExcluded = function(dvi, maxIncomp = 2, pairings = NULL, ignoreSex = FALSE, verbose = TRUE) {
  
  if(verbose) {
    cat("Finding exclusions\n")
    cat("Max incompatible markers =", maxIncomp, "\n")
  }
  
  # Ensure proper dviData object
  dvi = consolidateDVI(dvi)
  
  # Pairings
  if(!missing(pairings) || !missing(ignoreSex))
    pairings = pairings %||% generatePairings(dvi, ignoreSex = ignoreSex)
  else
    pairings = dvi$pairings %||% generatePairings(dvi, ignoreSex = ignoreSex)
  
  am = dvi$am
  missing = dvi$missing
  vics = names(dvi$pm)
  comp = getFamily(dvi, missing)
  
  # Main calculation: Exclusion matrix
  mat = exclusionMatrix(dvi, pairings = pairings)

  # Utilities for analysing rows/cols of `mat`
  minEx = function(v)
    if(all(is.na(v))) Inf else min(v, na.rm = TRUE)
  
  commnt = function(v) 
    ifelse(v == Inf, "No compatible pairings", sprintf("%s+ inconsistencies", v))
  
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
                           Conclusion = "Excluded", 
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
                           Conclusion = "Excluded",
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
  
  keepPairs = apply(mat[keepVics, keepMissing, drop = FALSE], 1, function(rw)
      c("*", keepMissing[!is.na(rw) & rw <= maxIncomp]), simplify = FALSE)
  
  dviRed$pairings = keepPairs

  nRemov = sum(lengths(dvi$pairings)) - sum(lengths(keepPairs))
  if(verbose)
    cat(sprintf("Pairings excluded in total: %d\n", nRemov))
  

  # Return list -------------------------------------------------------------

  excluded = list(sample = excludedVics, missing = excludedMissing, fam = excludedFams)
  summary = list(PM = summaryPM, AM = summaryAM)
  
  list(exclusionMatrix = mat, excluded = excluded, dviReduced = dviRed, 
       summary = summary)
}



#' @rdname findExcluded
#' @importFrom forrel findExclusions
#' @export
exclusionMatrix = function(dvi, pairings = NULL, ignoreSex = FALSE) {
  
  # Ensure proper dviData object
  dvi = consolidateDVI(dvi)
  
  am = dvi$am
  missing = dvi$missing
  vics = names(dvi$pm)
  
  # Pairings
  if(!missing(pairings) || !missing(ignoreSex))
    pairings = pairings %||% generatePairings(dvi, ignoreSex = ignoreSex)
  else
    pairings = dvi$pairings %||% generatePairings(dvi, ignoreSex = ignoreSex)
  
  # AM components (for use in output)
  comp = getFamily(dvi, missing)
  
  # Initialise matrix
  mat = matrix(NA_integer_, nrow = length(vics), ncol = length(missing), 
               dimnames = list(vics, missing))
  
  # Loop through each pair of victim vs missing
  for(vic in vics) {
    vicdata = dvi$pm[[vic]]
    compatMiss = setdiff(pairings[[vic]], "*")
    for(mis in compatMiss) {
      ref = am[[comp[mis]]]
      mat[vic, mis] = length(findExclusions(ref, id = mis, candidate = vicdata))
    }
  }
  
  mat
}

