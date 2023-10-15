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
  
  # Nonidentifiable ---------------------------------------------------------

  if(verbose)
    cat("Nonidentifiable missing persons" |> dashpad())
  
  non = findNonidentifiable(dvi)
  dvi = non$dviReduced
  summ = non$summary
  if(verbose)
    if(!is.null(summ)) print(summ) else cat("None\n")
  
  iter = 0
  summaries = list()
  
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
      summ = excl$summary$AM
      summaries = c(summaries, list(summ))
      nRemov = sum(excl$exclusionMatrix > maxIncomp, na.rm = TRUE)
      if(verbose)
        if(!is.null(summ)) print(summ) else cat(sprintf("Removed %d candidate pairings\n", nRemov))
    }
    
    if(verbose)
      cat("Undisputed, iteration" |> paste(iter) |> dashpad())
    

    # Undisputed --------------------------------------------------------------
    
    und = findUndisputed(dvi, threshold = threshold, limit = limit, 
                         check = FALSE, verbose = debug)
    if(dviEqual(und$dviReduced, dvi)) {
      if(verbose) cat("No change; breaking loop\n")
      break
    }
    else {
      dvi = und$dviReduced
      summ = und$summary
      summaries = c(summaries, list(summ))
      if(verbose)
        print(summ)
    }
  }
  
  # AM-driven ---------------------------------------------------------------

  if(verbose) {
    cat("AM-driven analysis" |> dashpad())
    nam = length(dvi$am)
    if(nam == 0) cat("0 remaining families\n")
    else if(nam == 1) cat("1 remaining family:", names(dvi$am), "\n\n")
    else cat(sprintf("%d remaining families: %s\n\n", nam, toString(names(dvi$am))))
  }
  
  amd = amDrivenDVI(dvi, threshold = threshold, threshold2 = threshold2, 
                    verbose = debug)
  
  dvi = amd$dviReduced
  summ = amd$summary
  summaries = c(summaries, list(summ))
  if(verbose && !is.null(summ))
    print(summ)
  

  # Remaining: Inconclusive -------------------------------------------------

  if(verbose)
    cat("Remaining MPs" |> dashpad())
    
  if(length(miss <- dvi$missing)) {
    summ = data.frame(Family = getFamily(dvi, miss),
                      Missing = miss,
                      Conclusion = "Inconclusive", 
                      row.names = NULL)
    summaries = c(summaries, list(summ))
    if(verbose)
      print(summ)
  }
  else {
    if(verbose) cat("None\n")
  }
  
  # Return final summary ----------------------------------------------------
    
  combineSummaries(summaries, orderBy = c("Family", "Missing"), dvi = origdvi)
}



dashpad = function(x, width = 50) {
  y = paste0("\n", strrep("-", 6), " ", x, " ")
  paste0(y, strrep("-", width - nchar(y)), "\n\n")
}