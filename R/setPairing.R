#' Set identifications manually
#'
#' Manually set one or several identifications in a DVI dataset. Typically,
#' these are obtained by external means, e.g., fingerprints, dental records etc.
#'
#' The command `setPairing(dvi, c("V" = "M"))` does the following:
#'
#' * Transfer the data of victim "V" to the individual "M" in the appropriate
#' reference family
#' * Remove "M" from the list of missing persons
#' * Remove "V" from the list of victim samples
#' * Update the list of pairings
#'
#' @param dvi A DVI dataset.
#' @param match A named vector of the format c(vic1 = miss2, vic2 = miss2, ...).
#' @param victim A vector of victim sample names. If NULL, defaulting to
#'   `names(match)`.
#' @param missing A vector of missing person names, of the same length as
#'   `victim`. If NULL, defaulting to `as.character(match)`.
#' @param Conclusion A character passed on to the `Conclusion` column of the
#'   output summary.
#' @param Comment A character passed on to the `Comment` column of the output
#'   summary.
#' @param verbose A logical, by default TRUE.
#'
#' @return A list with the following entries:
#'
#'   * `dviReduced`: The new `dviData` object, as described in Details
#'
#'   * `summary`: A data frame summarising the identifications
#'
#' @examples
#' x = setPairing(example2, match = c("V3" = "M2"))
#' x$dviReduced
#' x$summary
#'
#' # Alternative syntax, using `victim` and `missing`
#' y = setPairing(planecrash, victim = c("V4", "V5"), missing = c("M4", "M5"),
#'            Conclusion = "External evidence", Comment = "Dental")
#' y$dviReduced
#' y$summary
#' 
#' @export
setPairing = function(dvi, match = NULL, victim = NULL, missing = NULL, 
                      Conclusion = "Provided", Comment = "", verbose = TRUE) {
  dvi = consolidateDVI(dvi)
  
  if(!is.null(match)) {
    nms = names(match)
    if(is.null(nms))
      stop2("`match` must be a named vector: c(vic1 = miss1, ...)")
    if(!is.null(victim))
      stop2("When `match` is given, `victim` must be NULL, not: ", victim)
    if(!is.null(missing))
      stop2("When `match` is given, `missing` must be NULL, not: ", missing)
    
    victim = nms
    missing = as.character(match)
  }

  N = length(victim)
  if(N != length(missing))
    stop2(sprintf("Number of `victims` (%d) differs from the number of missing persons (%d)",
                  N, length(missing)))
  if(N == 0) {
    if(verbose) cat("No identifications indicated\n")
    return(list(dviReduced = dvi, summary = NULL))
  }
  
  if(length(badvic <- setdiff(victim, names(dvi$pm))))
    stop2("Unknown victim ID: ", badvic)
  if(length(badmis <- setdiff(missing, dvi$missing)))
    stop2("Unknown missing person ID: ", badmis)
  
  if(!length(Conclusion) %in% c(1, N))
    stop2(sprintf("`Conclusion` must have length either 1 or %d, not %d", N, length(Conclusion)))
  if(!length(Comment) %in% c(1, N))
    stop2(sprintf("`Comment` must have length either 1 or %d, not %d", N, length(Comment)))
  
  # AM component of each missing
  comp = getFamily(dvi, missing)
  
  if(verbose) {
    cat("Fixing the following identifications:\n")
    cat(sprintf(" %s = %s (Family %s)", victim, missing, comp), sep = "\n")
    cat("\n")
  }
  
  # Form reduced DVI Remove victims and missing
  dviRed = subsetDVI(dvi, pm = setdiff(names(dvi$pm), victim), 
                     missing = setdiff(dvi$missing, missing), verbose = verbose)
  
  # Transfer marker data
  dviRed$am = transferMarkers(from = dvi$pm[victim], to = dviRed$am, 
                              idsFrom = victim, idsTo = missing, erase = FALSE)
  
  # Build summary report
  summary = data.frame(Family = comp, Missing = missing, Sample = victim, 
                       Conclusion = Conclusion, Comment = Comment,
                       row.names = NULL)
  # Return
  list(dviReduced = dviRed,
       summary = summary)
}


#' Exclude pairings
#'
#' Disallow certain pairings by removing them from the list `dvi$pairings` of
#' candidate pairings for a given victim sample.
#' 
#' @param dvi A `dviData` object.
#' @param victim The name of a single victim sample.
#' @param missing The name(s) of one or several missing individuals.
#'
#' @return A `dviData` object.
#'
#' @examples
#' # Disallow V1 = M1 in the `example2` dataset:
#' ex = excludePairing(example2, victim = "V1", missing = "M1")
#' jointDVI(ex, verbose = FALSE)
#' 
#' # Compare with original
#' jointDVI(example2, verbose = FALSE)
#' 
#' # The only difference is in the `pairings` slot:
#' ex$pairings
#' example2$pairings
#' 
#' @export
excludePairing = function(dvi, victim, missing) {
  if(is.null(dvi$pairings))
    stop2("`dvi` object has no pairings attached")
  if(length(victim) != 1)
    stop2("Pairings can be excluded for only 1 victim at a time: ", victim)
  
  dvi$pairings[[victim]] = setdiff(dvi$pairings[[victim]], missing)
  dvi
}