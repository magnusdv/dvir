#' AM driven DVI
#'
#' Work-in-progress. Note: This function assumes that undisputed identifications
#' have been removed.
#'
#' @param dvi A `dviData` object.
#' @param fams A character; the names of families to consider. By default, all
#'   families. Special keywords: "simple" (all families with exactly 1 missing)
#'   and "nonsimple" (all families with > 1 missing).
#' @param threshold LR threshold for 'certain' match.
#' @param threshold2 LR threshold for 'probable' match (in *simple* families).
#' @param removeIf A character; identifications whose conclusion statement
#'   matches one of these (through `grepl`) are removed in the reduced dataset.
#'   By default: "undisputed", "probable" and "no match".
#' @param verbose A logical.
#'
#' @return A list of `dviReduced` and `summary`.
#'
#' @examples
#' w = amDrivenDVI(example2)
#' w$summary
#' w$dviReduced
#' 
#' # Bigger example: Undisputed first
#' u = findUndisputed(planecrash)
#' u$summary
#'
#' # AM-driven analysis of the remaining
#' amDrivenDVI(u$dviReduced, threshold2 = 500)
#'
#' @export
amDrivenDVI = function(dvi, fams = NULL, threshold = 1e4, threshold2 = 10^3, 
                       removeIf = c("undisputed", "probable", "no match"),
                       verbose = TRUE) {
  dvi = consolidateDVI(dvi)
  
  if(verbose)
    cat("AM-driven analysis\n")
  
  if(is.null(fams))
    famnames = names(dvi$am)
  else if(identical(fams, "simple"))
    famnames = getSimpleFams(dvi)
  else if(identical(fams, "nonsimple"))
    famnames = setdiff(names(dvi$am), getSimpleFams(dvi))
  else
    famnames = fams

    
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
  
  # Subset of summary with fixed conclusion
  patt = paste(removeIf, collapse = "|")
  s = summary[grepl(patt, tolower(summary$Conclusion)), , drop = FALSE]
  
  if(is.null(summary) || nrow(s) == 0) {
    if(verbose)
      cat("No reduction of the dataset\n")
    return(list(dviReduced = dvi, summary = summary))
  }
  
  remainMissing = setdiff(dvi$missing, unlist(strsplit(s$Missing, split = ",")))
  remainVics = setdiff(names(dvi$pm), unlist(strsplit(s$Sample, split = ",")))
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
    comment = paste("Also:", paste(sprintf("%s (LR=%.2g)", names(also), also), collapse = ", "))
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


.jointFamDVI = function(dvi1, threshold, check = FALSE, verbose = FALSE) {
  fam = names(dvi1$am)
  missing = dvi1$missing
  
  j = jointDVI(dvi1, undisputed = FALSE, check = check, verbose = verbose)
  nrw = nrow(j)
  
  # LR column
  lrs = j$LR
  
  # Jointy undisputed: LR_1 >= thresh AND LR_1:2 >= thresh
  if(lrs[1] >= threshold && (nrw == 1 || lrs[1]/lrs[2] >= threshold)) {
    # Compactify joint data frame
    res0 = compactJointRes(j[1, ])
    
    # Remove columns with '*' (should not be reported)
    goodcols = apply(res0, 2, function(cc) all(cc %in% missing))
    res = res0[, goodcols, drop = FALSE]
    
    vics = names(res)
    miss = as.character(res)

    # Character with identified pairs
    prs = sprintf("%s=%s", miss, vics)
    
    summary = data.frame(Family = fam, Missing = miss, Sample = vics, LR = lrs[1],
               Conclusion = "Jointly undisputed",
               Comment = paste("Joint with:", sapply(seq_along(prs), function(i) toString(prs[-i]))))
    return(summary)
  }
  
  if(nrw < 3)
    return(NULL)

  if(lrs[1] == lrs[2] && lrs[1] >= threshold/2 && lrs[1]/lrs[3] >= threshold) {
    # Symmetric pair of solutions (e.g. indistinguishable siblings)
    
    # Compactify joint data frame
    res0 = compactJointRes(j[1:2, ])
    
    # Remove columns with '*' (should not be reported)
    goodcols = apply(res0, 2, function(cc) all(cc %in% missing))
    res = res0[, goodcols, drop = FALSE]
    
    if(!setequal(res[1,], res[2,])) {
      message("Warning: This type of symmetry is currently only partially reported:")
      print(res0)
      
      eq = res[1,] == res[2, ]
      vics = names(res)[eq]
      miss = as.character(res[1, eq])
      conc = "Jointly undisputed"
      cmt = paste("Joint with", sapply(seq_along(miss), function(i) toString(miss[-i])))
    }
    else {
      vics = names(res) |> paste(collapse = "/")
      miss = as.character(res[1,])
      conc = "Symmetric undisputed"
      cmt = paste("Symmetric with", sapply(seq_along(miss), function(i) toString(miss[-i])))
    }
    
    summary = data.frame(Family = fam, Missing = miss, Sample = vics, 
                         LR = lrs[1] * 2, Conclusion = conc, Comment = cmt)
  }
  else {
    summary = NULL
  }
  
  summary
}