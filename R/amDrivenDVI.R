
#' AM driven DVI
#'
#' Work-in-progress. Note: This function assumes that undisputed identifications
#' have been removed.
#' 
#' @param dvi A `dviData` object.
#' @param simple A logical (only TRUE for now)
#' @param thresholdProbable LR treshold
#' @param LRmatrix NULL
#' @param ignoreSex FALSE
#' @param verbose TRUE
#'
#' @return A list of `dviReduced` and `summary`.
#' 
#' @examples
#' amDrivenDVI(example2)
#' 
#' # Bigger example: Undisputed first
#' u = findUndisputed(planecrash)
#' u$summary
#' amDrivenDVI(u$dviReduced)
#' 
#' @export
amDrivenDVI = function(dvi, simple = TRUE, thresholdProbable = 10^3, 
                       LRmatrix = NULL, ignoreSex = FALSE, verbose = TRUE) {
  dvi = consolidateDVI(dvi)
  
  # Enforce family names
  if(is.null(names(dvi$am)))
    dvi = relabelDVI(dvi, familyPrefix = "")
  
  if(!simple)
    stop2("Only simple families are supperted for now")
  
  simple = getSimpleFams(dvi)
  if(!length(simple))
    stop2("No simple reference families in the dataset\n")
  if(verbose)
    cat("Identified simple families:", toString(simple), "\n")
  
  if(is.null(LRmatrix)) {
    .dviRed = subsetDVI(dvi, am = simple, verbose = FALSE)
    LRmatrix = pairwiseLR(.dviRed, verbose = FALSE)$LRmatrix
  }
  
  resList = lapply(simple, function(s) {
    dvi1 = subsetDVI(dvi, am = s, verbose = FALSE)
    if(verbose)
      cat(sprintf("Analysing family '%s' (missing person '%s')", s, dvi1$missing), "\n")
    r = .simpleFamDVI(dvi1, threshold = thresholdProbable, LRmatrix = LRmatrix)
    r
  })
  
  summary = do.call(rbind, resList)
  
  nonsimple = setdiff(names(dvi$am), simple)
  nam = length(nonsimple)
  if(nam > 0)
    dviRed = subsetDVI(dvi, am = nonsimple, verbose = FALSE)
  else 
    dviRed = NULL
  
  if(verbose) 
    cat("Reduced dataset:", if(!nam) "Empty!" else paste(nam, if(nam==1) "family" else "families"), "\n")

  list(dviReduced = dviRed, summary = summary)
}

.simpleFamDVI = function(dviSimple, threshold, LRmatrix = NULL, ignoreSex = FALSE) {
  
  miss = dviSimple$missing
  if(length(miss) != 1)
    stop2("`.simpleFamDVI()` expects a dataset with a single missing person: ", miss)
  
  am = dviSimple$am
  if(length(am) != 1 || is.null(names(am))) {
    print(dviSimple)
    stop2("`.simpleFamDVI()` expects a `dviData` object with a single, named, family")
  }
  
  if(is.null(LRmatrix))
    LRmatrix = pairwiseLR(dviSimple, ignoreSex = ignoreSex, 
                          check = FALSE, verbose = FALSE)$LRmatrix
  
  # Row of LRs for the missing person
  lrs = LRmatrix[, miss]
  names(lrs) = rownames(LRmatrix)
  
  # Max LR and the corresponding victim (first of, if several)
  maxLR = max(lrs)
  bestMatch = names(lrs)[which.max(lrs)]
  
  # All LRs exceeding threshold
  top = lrs[lrs > threshold]
  comment = ""
  
  if(length(top) > 1) {
    concl = "Disputed"
    comment = paste(sprintf("%s (LR=%.2g)", names(top), top), collapse = ", ")
  }
  else if(length(top) == 1) {
    concl = "Probable"
    lrs2 = lrs[names(lrs) != bestMatch]
    if(length(lrs2)) {
      runnerup = names(lrs2)[which.max(lrs2)]
      comment = sprintf("Runner-up: %s (LR=%.2g)", runnerup, max(lrs2))
    }
    else 
      comment = "Runner-up: -"
  }
  else {
    concl = "No match"
    comment = sprintf("Best: %s (LR=%.2g)", bestMatch, maxLR)
    # Blanks in summary
    bestMatch = NA_character_
    maxLR = NA_real_
    maxLR = NA
  }
  
  data.frame(Family = names(dviSimple$am), Missing = miss, Sample = bestMatch, 
             LR = maxLR, Conclusion = concl, Comment = comment)
}