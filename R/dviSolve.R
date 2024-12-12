#' A complete pipeline for solving a DVI case
#'
#' This wraps several other functions into a complete pipeline for solving a DVI
#' case.
#'
#' @param dvi A `dviData` object.
#' @param threshold LR threshold for 'significant' match.
#' @param threshold2 LR threshold for 'probable' match. By default set to
#'   `threshold/10`.
#' @param maxIncomp	An integer passed onto [findExcluded()]. A pairing is
#'   excluded if the number of incompatible markers exceeds this.
#' @param ignoreSex A logical, by default FALSE.
#' @param limit	A number passed onto [findUndisputed()]; only pairwise LR values
#'   above this are considered.
#' @param detailedOutput A logical, by default FALSE. See Details.
#' @param verbose,debug Logicals.
#'
#' @return A list of data frames `AM` and `PM`. If `detailedOutput` is TRUE, the
#'   LR matrix and exclusion matrix from the first iteration are also included.
#'
#' @examples
#' dviSolve(example2)
#' dviSolve(example2, threshold = 5, detailedOutput = TRUE, verbose = FALSE)
#'
#' @export
dviSolve = function(dvi, threshold = 1e4, threshold2 = max(1, threshold/10), 
                    maxIncomp = 2, ignoreSex = FALSE, limit = 0, 
                    detailedOutput = FALSE, verbose = TRUE, debug = FALSE) {

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
  summariesAM = c(summariesAM, list(non$summary))
  if(verbose)
    printSummary(non$summary, nulltext = "None")
  
  iter = 0

  # Loop until no change: excl, undisp, excl, ...
  while(TRUE) {
    
    iter = iter + 1
    
    # Exclusions --------------------------------------------------------------
    
    if(verbose)
      cat("Exclusions, iteration" |> paste(iter) |> dashpad())
    
    excl = findExcluded(dvi, maxIncomp = maxIncomp, verbose = debug)

    # Store first exclusion matrix
    if(iter == 1)
      EXCLmat1 = excl$exclusionMatrix
    
    # Break loop if no change - but not in first iter!
    if(dviEqual(excl$dviReduced, dvi)) {
      if(verbose) cat("No exclusions", if(iter > 1) "; breaking loop", "\n", sep = "")
      if(iter > 1)
        break
    }
    else {
      s = excl$summary
      summariesAM = c(summariesAM, list(s$AM))
      summariesPM = c(summariesPM, list(s$PM))
      if(verbose) {
        removedPrs = removedPairings(excl$dviReduced, dvi)
        printSummary(s, addNewline = removedPrs > 0)
        
        if(removedPrs)
          cat(sprintf("Excluded %s\n", pluralise("individual pairing", removedPrs)))
      }
      dvi = excl$dviReduced
    }
    
    # Undisputed --------------------------------------------------------------
    
    if(verbose)
      cat("Undisputed, iteration" |> paste(iter) |> dashpad())
    
    und = findUndisputed(dvi, threshold = threshold, limit = limit, 
                         keepLRmatrs = (iter == 1), verbose = debug)

    # Store first LR matrix
    if(iter == 1) {
      LRmat1 = und$LRmatrix[[1]]
      LRmat = und$LRmatrix[[length(und$LRmatrix)]]
    }
    else 
      LRmat = und$LRmatrix
    
    if(dviEqual(und$dviReduced, dvi)) {
      if(verbose) cat("No change; breaking loop\n")
      break
    }
    else {
      s = und$summary
      summariesAM = c(summariesAM, list(s))
      summariesPM = c(summariesPM, list(s))
      if(verbose) {
        removedPrs = removedPairings(und$dviReduced, dvi)
        printSummary(s, addNewline = removedPrs > 0)
        if(removedPrs) {
          a = pluralise("pairing", removedPrs)
          if(limit == 0) 
            cat(sprintf("Removed %s with LR = 0\n", a))
          else 
            cat(sprintf("Removed %s with LR < %g (`limit`)\n", a, limit))
        }
      }
      dvi = und$dviReduced
    }
  }

  # AM-driven: Simple-----------------------------------------------------------
  
  nam = length(dvi$am)
  nMiss = nMissFam(dvi)
  simpleFams = names(nMiss[nMiss == 1])
  nsimp = length(simpleFams)
  
  if(verbose) {
    cat("AM-driven analysis: Simple families" |> dashpad())
    if(nsimp == 0)
      cat("0 simple families remaining\n")
    else
      cat(sprintf("%d simple famil%s remaining: %s\n\n", nsimp, ifelse(nsimp == 1, "y", "ies"), toString(simpleFams)))
  }
  
  for(fam in simpleFams) {
    dvi1 = subsetDVI(dvi, am = fam, verbose = FALSE)
    if(verbose)
      cat(sprintf("Family %s: %s", fam, dvi1$missing))
    s = .simpleFamDVI(dvi1, threshold = threshold2, LRmatrix = LRmat)
    summariesAM = c(summariesAM, list(s))
    if(verbose)
      cat(" -->", s$Conclusion, "\n")
  }
  
  # AM-driven: Complex-----------------------------------------------------------
  
  complexFams = names(nMiss[nMiss > 1])
  ncomp = length(complexFams)
  
  if(verbose) {
    cat("AM-driven analysis: Complex families" |> dashpad())
    if(ncomp == 0)
      cat("0 complex families remaining\n")
    else
      cat(sprintf("%d complex famil%s remaining: %s\n\n", ncomp, 
                  ifelse(ncomp == 1, "y", "ies"), toString(complexFams)))
  }

  # Loop through families; remove identified victims in each iteration
  for(fam in complexFams) {
    dvi1 = subsetDVI(dvi, am = fam, verbose = FALSE)
    if(verbose)
      cat(sprintf("Family %s: %s", fam, toString(dvi1$missing)))
    s = .complexFamDVI(dvi1, threshold = threshold, LRmatrix = LRmat, verbose = FALSE)
    summariesAM = c(summariesAM, list(s$AM))
    summariesPM = c(summariesPM, list(s$PM))

    # Reduce main DVI dataset
    dvi = subsetDVI(dvi, 
                    am = .mysetdiff(names(dvi$am), fam),
                    pm = .mysetdiff(names(dvi$pm), s$PM$Sample), 
                    removeUnpairedPM = FALSE, verbose = FALSE)
    if(verbose) {
      concs = s$AM$Conclusion
      if(all(concs == concs[1])) concs = concs[1]
      cat(" -->", toString(concs), "\n")

    }
  }
  
  # PM driven analysis----------------------------------------------------------
  
  if(verbose)
    cat("Remaining victim samples" |> dashpad())
  
  vics = names(dvi$pm)
  nv = length(vics)
  if(verbose) {
    if(nv == 0) cat("No remaining victims\n")
    else cat(sprintf("%d remaining victim%s: %s\n", nv, if(nv == 1) "" else "s", toString(vics)))
  }

  if(nv > 0) {
    s = .pmDrivenDVI(dvi, threshold2, LRmatrix = LRmat, origdvi = origdvi)
    summariesPM = c(summariesPM, list(s))
  }

  # Collect results
  resultAM = formatSummary(summariesAM, orientation = "AM", dvi = origdvi)
  resultPM = formatSummary(summariesPM, orientation = "PM", dvi = origdvi)

  res = list(AM = resultAM, PM = resultPM)
  if(detailedOutput) {
    res$LRmatrix = LRmat1
    res$exclusionMatrix = EXCLmat1
  }
  
  res
}



dashpad = function(x, width = 50) {
  y = paste0("\n", strrep("-", 6), " ", x, " ")
  paste0(y, strrep("-", width - nchar(y)), "\n\n")
}

printSummary = function(s, nulltext = NULL, addNewline = FALSE) {
  if(is.null(s) || (length(s) == 2 && is.null(s$AM) && is.null(s$PM))) {
    if(!is.null(nulltext))
      cat(nulltext, "\n")
    return()
  }
  if(is.data.frame(s))
    print(s)
  if(!is.null(s$AM)) {
    cat("$AM\n"); print(s$AM)
  }
  if(!is.null(s$PM)) {
    cat("$PM\n"); print(s$PM)
  }
  if(addNewline)
    cat("\n")
}