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
  
  # AM components
  comp = getComponent(am, missing, checkUnique = TRUE, errorIfUnknown = TRUE)
  
  # Loop through each pair of victim vs missing
  for(i in 1:npm) for(j in 1:nmiss) {
    vic = pm[[i]]
    ref = am[[comp[j]]]
    mat[i, j] = length(findExclusions(ref, id = missing[j], candidate = vic))
  }
  
  mat
}


#' Excluded pairings
#'
#' @param dvi A [dviData()] object.
#' @param pairings A list of possible pairings for each victim. If NULL, all
#'   sex-consistent pairings are used, unless if `ignoreSex` is set to TRUE (in
#'   which case all pairings are considered).
#' @param ignoreSex A logical, by default FALSE.
#' @param maxIncomp An integer. A pairing is excluded if the number of
#'   incompatible markers exceeed this.
#' @param verbose A logical. Default: TRUE.
#'
#' @return A list with the following entries:
#'
#'   * `exclusionMatrix`: A matrix showing the number of inconsistencies for
#'   each pair (or NA if the pairing was not considered)
#'
#'   * `pmNomatch`: A vector with the names of victims (PM singletons) who are
#'   excluded against all missing persons.
#'
#'   * `missNomatch`: A vector with the missing persons who are excluded against
#'   all victims.
#'
#'   * `dviReduced`: A reduced version of `dvi`, where excluded victims/missing
#'   persons are removed.
#'   
#'   * `pairings` A list of non-excluded pairings.
#'
#' @examples
#' findExcluded(icmp)
#'
#' @export
findExcluded = function(dvi, pairings = NULL, ignoreSex = FALSE, maxIncomp = 2, verbose = TRUE) {
  
  if(verbose) {
    message("Finding exclusions")
    message("Max incompatible markers = ", maxIncomp)
  }
  
  # Ensure proper dviData object
  dvi = consolidateDVI(dvi)
  
  if(is.null(pairings))
    pairings = generatePairings(dvi, ignoreSex = ignoreSex)
  
  pm = dvi$pm
  am = dvi$am
  missing = dvi$missing
  vics = names(pm)
  

  # Exclusion matrix --------------------------------------------------------

  # Initialise matrix
  mat = matrix(NA_integer_, nrow = length(pm), ncol = length(missing), 
               dimnames = list(vics, missing))
  
  # AM components
  comp = getComponent(am, missing, checkUnique = TRUE, errorIfUnknown = TRUE)
  names(comp) = missing
  
  # Loop through each pair of victim vs missing
  for(vic in vics) {
    vicdata = pm[[vic]]
    compatMiss = setdiff(pairings[[vic]], "*")
    for(mis in compatMiss) {
      ref = am[[comp[mis]]]
      mat[vic, mis] = length(findExclusions(ref, id = mis, candidate = vicdata))
    }
  }

  # Victims excluded against all --------------------------------------------
  
  pmNomatch = vics[apply(mat, 1, function(a) all(is.na(a) | a > maxIncomp))]
  
  if(verbose) {
    message("\nPM samples excluded against all missing:", if(!length(pmNomatch)) " None")
    for(m in pmNomatch)
      message(sprintf(" %s (minimum %s inconsistencies)", m, min(mat[m, ], na.rm = TRUE)))
  }
  
  # Missing excluded against all --------------------------------------------
  
  missNomatch = missing[apply(mat, 2, function(a) all(is.na(a) | a > maxIncomp))]
  
  if(verbose) {
    message("\nMissing persons excluded against all PM samples:", if(!length(missNomatch)) " None")
    for(m in missNomatch)
      message(sprintf(" %s (minimum %s inconsistencies)", m, min(mat[, m], na.rm = TRUE)))
  }

  # Reduced DVI problem -----------------------------------------------------
  
  keepVics = setdiff(vics, pmNomatch)
  keepMissing = setdiff(missing, missNomatch)
  message("")
  if(length(pmNomatch) + length(missNomatch))
    dviRed = subsetDVI(dvi, pm = keepVics, missing = keepMissing, verbose = verbose)
  else {
    if(verbose) 
      message("No reduction of the DVI dataset")
    dviRed = dvi
  }
  
  # Updated pairings ----------------------------------------------------
  
  nRemov = sum(mat > maxIncomp, na.rm = TRUE)
  if(nRemov > 0) {
    if(verbose)
      message(sprintf("\nIn total %d pairings excluded", nRemov))
    
    keepPairs = apply(mat[keepVics, , drop = FALSE], 1, function(rw)
      c("*", missing[!is.na(rw) & rw <= maxIncomp]), simplify = FALSE)
  }
  else {
    if(verbose) message("No pairings excluded")
    keepPairs = pairings
  }
    
  # Return list -------------------------------------------------------------

  list(exclusionMatrix = mat, pmNomatch = pmNomatch, missNomatch = missNomatch, 
       dviReduced = dviRed, pairings = keepPairs)
}
