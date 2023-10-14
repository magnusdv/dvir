
#' AM driven DVI
#'
#' Work-in-progress. Note: This function assumes that undisputed identifications
#' have been removed.
#' 
#' @param dvi A `dviData` object.
#' @param threshold LR threshold for 'certain' match 
#' @param threshold2 LR threshold for 'probable' match
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
#' 
#' # AM-driven analysis of the remaining
#' amDrivenDVI(u$dviReduced, threshold2 = 500)
#' 
#' @export
amDrivenDVI = function(dvi, threshold = 1e4, threshold2 = 10^3, verbose = TRUE) {
  dvi = consolidateDVI(dvi)
  
  if(verbose)
    cat("AM-driven analysis\n")
  
  # Enforce family names   # TODO: Move to consolidateDVI()?
  if(is.null(names(dvi$am)))
    dvi = relabelDVI(dvi, familyPrefix = "")
  
  famnames = names(dvi$am)
  
  # Number of missing in each fam
  comp = getFamily(dvi, ids = dvi$missing)
  nMiss = sapply(famnames, function(fam) sum(comp == fam))
  names(nMiss) = famnames
  
  # If any simple families, compute LR matrix once
  
  if(any(nMiss == 1)) {
    simple = names(which(nMiss == 1))
    .dviRed = subsetDVI(dvi, am = simple, verbose = FALSE)
    LRmatrix = pairwiseLR(.dviRed, verbose = FALSE)$LRmatrix
  }
  
  resList = lapply(famnames, function(fam) {
    dvi1 = subsetDVI(dvi, am = fam, verbose = FALSE)
    if(verbose)
      cat(sprintf(" Family %s; %d missing (%s)", fam, length(dvi1$missing), toString(dvi1$missing)))
    if(nMiss[fam] == 1)
      r = .simpleFamDVI(dvi1, threshold = threshold2, LRmatrix = LRmatrix)
    else
      r = .jointFamDVI(dvi1, threshold = threshold)
    if(verbose)
      cat(" -->", if(is.null(r)) "Inconclusive" else r$Conclusion[1], "\n")

    r
  })
  
  summary = do.call(rbind, resList)
  
  remainMissing = setdiff(dvi$missing, summary$Missing)
  remainVics = setdiff(names(dvi$pm), summary$Sample)
  if(length(remainMissing) || length(remainVics))
    dviRed = subsetDVI(dvi, pm = remainVics, missing = remainMissing, verbose = FALSE)
  else 
    dviRed = NULL
  
  if(verbose)  {
    nam = length(dviRed$am)
    cat("Reduced dataset:", if(nam==0) "Empty!\n" else paste(nam, if(nam==1) "family" else "families", "remaining\n"))
  }
  
  list(dviReduced = dviRed, summary = summary)
}

.simpleFamDVI = function(dvi1, threshold, LRmatrix = NULL) {
  
  miss = dvi1$missing
  if(length(miss) != 1)
    stop2("`.simpleFamDVI()` expects a dataset with a single missing person: ", miss)
  
  am = dvi1$am
  if(length(am) != 1 || is.null(names(am))) {
    print(dvi1)
    stop2("`.simpleFamDVI()` expects a `dviData` object with a single, named, family")
  }
  
  if(is.null(LRmatrix))
    LRmatrix = pairwiseLR(dvi1, check = FALSE, verbose = FALSE)$LRmatrix
  
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
    also = top[names(top) != bestMatch]
    comment = paste(sprintf("%s (LR=%.2g)", names(also), also), collapse = ", ")
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
  
  data.frame(Family = names(dvi1$am), Missing = miss, Sample = bestMatch, 
             LR = maxLR, Conclusion = concl, Comment = comment)
}


.jointFamDVI = function(dvi1, threshold, ...) {
  fam = names(dvi1$am)
  
  j = jointDVI(dvi1, undisputed = FALSE, ..., verbose = FALSE)
  
  # Require both LR_1 >= threshold AND LR_1:2 >= threshold
  lrs = j$LR
  good = lrs[1] >= threshold && lrs[1]/lrs[2] >= threshold
  
  if(good) {
    # Compactify joint data frame
    res = compactJointRes(j[1, ])
    
    # Number of pairs in solution
    n = ncol(res) - 3 
    vics = names(res)[seq_len(n)]
    miss = as.character(res[1, seq_len(n)])

    # Character with identified pairs
    prs = sprintf("%s=%s", miss,vics)
    
    summary = data.frame(Family = fam, Missing = miss, Sample = vics, LR = lrs[1],
               Conclusion = "Jointly undisputed",
               Comment = paste("Joint with", sapply(seq_along(prs), function(i) toString(prs[-i]))))
  }
  else {
    summary = NULL
  }
  
  summary
}