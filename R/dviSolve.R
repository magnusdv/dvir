#' A complete pipeline for solving a DVI case
#'
#' This wraps several other functions into a complete pipeline for solving a DVI
#' case.
#' 
#' @param dvi A `dviData` object.
#' @param threshold LR threshold for 'significant' match.
#' @param threshold2 LR threshold for 'probable' match.
#' @param maxIncomp	An integer passed onto [findExcluded()]. A pairing is
#'   excluded if the number of incompatible markers exceeds this.
#' @param ignoreSex A logical, by default FALSE.
#' @param limit	A number passed onto [findUndisputed()]; only pairwise LR values
#'   above this are considered.
#' @param verbose,debug Logicals.
#'
#' @return A data frame.
#'
#' @examples
#' dviSolve(example2)
#' dviSolve(example2, threshold = 5, verbose = FALSE)
#' @export
dviSolve = function(dvi, threshold = 1e4, threshold2 = max(1, threshold/10), maxIncomp = 2, 
                    ignoreSex = FALSE, limit = 0, verbose = TRUE, debug = FALSE) {

  if(ignoreSex)
    dvi$pairings = generatePairings(dvi, ignoreSex = TRUE)
  
  # Check dataset -----------------------------------------------------------
  if(verbose)
    cat("Checking dataset" |> dashpad())
  
  origdvi = dvi = consolidateDVI(dvi)
  checkDVI(dvi, verbose = debug, ignoreSex = ignoreSex)
  if(verbose)
    cat("Ok\n")
  
  summariesAM = summariesPM = list()
  
  # Nonidentifiable ---------------------------------------------------------

  if(verbose)
    cat("Nonidentifiable missing persons" |> dashpad())
  
  non = findNonidentifiable(dvi)
  dvi = non$dviReduced
  summ = non$summary
  summariesAM = c(summariesAM, list(summ))
  if(verbose)
    if(!is.null(summ)) print(summ) else cat("None\n")
  
  iter = 0

  # Loop until no change: excl, undisp, excl, ...
  while(TRUE) {
    iter = iter + 1
    
    # Exclusions --------------------------------------------------------------
    
    if(verbose)
      cat("Exclusions, iteration" |> paste(iter) |> dashpad())
    
    excl = findExcluded(dvi, maxIncomp = maxIncomp, verbose = debug)
    if(dviEqual(excl$dviReduced, dvi)) {
      if(verbose) cat("No change; breaking loop\n")
      break
    }
    else {
      dvi = excl$dviReduced
      summ = excl$summary
      summariesAM = c(summariesAM, list(summ$AM))
      summariesPM = c(summariesPM, list(summ$PM))
      nRemov = sum(excl$exclusionMatrix > maxIncomp, na.rm = TRUE)
      if(verbose) {
        if(!is.null(summ)) {print(summ); cat("\n")}
        cat(sprintf("Removed %d candidate pairings\n", nRemov))
      }
    }
    
    if(verbose)
      cat("Undisputed, iteration" |> paste(iter) |> dashpad())
    

    # Undisputed --------------------------------------------------------------
    
    und = findUndisputed(dvi, threshold = threshold, limit = limit, verbose = debug)
    if(dviEqual(und$dviReduced, dvi)) {
      if(verbose) cat("No change; breaking loop\n")
      break
    }
    else {
      dvi = und$dviReduced
      summ = und$summary
      summariesAM = c(summariesAM, list(summ))
      summariesPM = c(summariesPM, list(summ))
      if(verbose)
        print(summ)
    }
  }

  # AM-driven ---------------------------------------------------------------

  nam = length(dvi$am)
  if(verbose) {
    cat("AM-driven analysis" |> dashpad())
    if(nam == 0) cat("0 remaining families\n")
    else if(nam == 1) cat("1 remaining family:", names(dvi$am), "\n\n")
    else cat(sprintf("%d remaining families: %s\n\n", nam, toString(names(dvi$am))))
  }
  
  if(nam >= 1) { 
    amd = amDrivenDVI(dvi, threshold = threshold, threshold2 = threshold2, 
                      verbose = debug)
    
    dvi = amd$dviReduced
    summariesAM = c(summariesAM, list(amd$summary$AM))
    summariesPM = c(summariesPM, list(amd$summary$PM))
    if(verbose && !is.null(amd$summary$AM))
      print(amd$summary$AM)
  }

  # Remaining: Inconclusive -------------------------------------------------

  if(verbose)
    cat("Remaining MPs" |> dashpad())
    
  if(length(miss <- dvi$missing)) {
    summAM = data.frame(Family = getFamily(dvi, miss),
                        Missing = miss,
                        Conclusion = "Inconclusive", 
                        row.names = NULL)
    summariesAM = c(summariesAM, list(summAM))
    if(verbose)
      print(summAM)
  }
  else if(verbose) 
    cat("None\n")
  
  
    if(verbose)
    cat("Remaining victim samples" |> dashpad())
    
  if(length(dvi$pm)) {
    summPM = data.frame(Sample = names(dvi$pm),
                        Conclusion = "Inconclusive", 
                        row.names = NULL)
    summariesPM = c(summariesPM, list(summPM))
    if(verbose)
      print(summPM)
  }
  else if (verbose)
    cat("None\n")
  
  # Return final summaries ----------------------------------------------------

  resultAM = formatSummary(summariesAM, centricity = "AM", dvi = origdvi)
  resultPM = formatSummary(summariesPM, centricity = "PM", dvi = origdvi)
  list(AM = resultAM, PM = resultPM)
}



dashpad = function(x, width = 50) {
  y = paste0("\n", strrep("-", 6), " ", x, " ")
  paste0(y, strrep("-", width - nchar(y)), "\n\n")
}
