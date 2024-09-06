#' AM driven DVI
#'
#' AM-driven identification, i.e., considering one AM family at a time. Simple
#' families (exactly 1 missing) are handled directly from the LR matrix, while
#' nonsimple families are analysed with [dviJoint()].
#'
#' Note: This function assumes that undisputed identifications have been
#' removed. Strange outputs may occur otherwise.
#'
#' @param dvi A `dviData` object.
#' @param fams A character; the names of families to consider. By default, all
#'   families. Special keywords: "simple" (all families with exactly 1 missing)
#'   and "nonsimple" (all families with > 1 missing).
#' @param threshold LR threshold for 'certain' match.
#' @param threshold2 LR threshold for 'probable' match (in *simple* families).
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
amDrivenDVI = function(dvi, fams = NULL, threshold = 1e4, threshold2 = max(1, threshold/10), 
                       verbose = TRUE) {
  dvi = consolidateDVI(dvi)
  
  if(verbose)
    cat("AM-driven analysis\n")
  
  if(!length(dvi$pm) || !length(dvi$missing))
    return(list(dviReduced = dvi, summary = NULL))
  
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
  
  # Compute LR matrix once and for all
  LRmat = pairwiseLR(dvi, check = FALSE, verbose = FALSE)$LRmatrix
  
  # List of summaries (to be combined below)
  summariesAM = summariesPM = NULL
  
  for(fam in famnames) {
    
    dvi1 = subsetDVI(dvi, am = fam, verbose = FALSE)
    if(verbose)
      cat(sprintf(" Family %s; %d missing (%s)", fam, length(dvi1$missing), 
                  toString(dvi1$missing)))
    
    if(nMiss[fam] == 1) {
      s0 = .simpleFamDVI(dvi1, threshold = threshold2, LRmatrix = LRmat)
      s = list(AM = s0, PM = s0)
    }
    else {
      s = .jointFamDVI(dvi1, threshold = threshold, LRmatrix = LRmat, verbose = FALSE)
    }
    
    summariesAM = c(summariesAM, list(s$AM))
    summariesPM = c(summariesPM, list(s$PM))
    if(verbose)
      cat(" -->", if(is.null(s$AM)) "Inconclusive" else toString(unique(s$AM$Conclusion)), "\n")
  }
  
  summaryAM = formatSummary(summariesAM, "AM")
  summaryPM = formatSummary(summariesPM, "PM")
  if(is.null(summaryAM)) {
    if(verbose)
      cat("No reduction of the dataset\n")
    return(list(dviReduced = dvi, summary = NULL))
  }
  
  remainMissing = setdiff(dvi$missing, summaryAM$Missing)
  remainVics = setdiff(names(dvi$pm), summaryPM$Sample)
  if(length(remainMissing) || length(remainVics))
    dviRed = subsetDVI(dvi, pm = remainVics, missing = remainMissing, verbose = FALSE)
  else 
    dviRed = NULL
  
  if(verbose)  {
    nam = length(dviRed$am)
    cat("Reduced dataset:", if(nam==0) "Empty!\n" else paste(nam, if(nam==1) "family" else "families", "remaining\n"))
  }
  
  list(dviReduced = dviRed, summary = list(AM = summaryAM, PM = summaryPM))
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
  
  # Summary
  data.frame(Family = names(dvi1$am), Missing = miss, Sample = bestMatch, 
             LR = maxLR, Conclusion = concl, Comment = comment)
}


.jointFamDVI = function(dvi1, threshold = 1e4, LRmatrix = NULL, verbose = FALSE, progress = verbose) {
  
  summaryAM = summaryPM = NULL
  
  # Joint table
  jres = dviJoint(dvi1, verbose = verbose, progress = progress)
  j = jres$joint %||% jres # TODO: clean up
  
  # Pairwise significant GLR
  pairGLR = pairwiseGLR(dvi1, jointTable = j, LRmatrix = LRmatrix, threshold = threshold)
  if(!is.null(s <- pairGLR$summary)) {
    summaryAM = summaryPM = s
    dvi1 = pairGLR$dviReduced
    j = j[, !names(j) %in% c(s$Sample, s$Missing), drop = FALSE]
  } 
  
  # Symmetric GLR
  symGLR = symmetricGLR(dvi1, jointTable = j, threshold = threshold, verbose = verbose)
  if(!is.null(s <- symGLR)) {
    summaryAM = formatSummary(list(summaryAM, s$AM), "AM")
    summaryPM = formatSummary(list(summaryPM, s$PM), "PM")
  }
  
  # Return summaries for this family
  list(AM = summaryAM, PM = summaryPM)
}

#' @importFrom verbalisr verbalise
symmetricGLR = function(dvi, jointTable = NULL, threshold = 1e4, verbose = FALSE) {
  fam = names(dvi$am)
  vics = names(dvi$pm)
  missing = dvi$missing
  j = jointTable %||% dviJoint(dvi, verbose = verbose)
  if(nrow(j) < 3)
    return(NULL)
  
  # Diffs between 1st and 2nd rows
  r1 = j[1, vics]; r2 = j[2, vics]
  diffs = r1 != r2
  if(sum(diffs) > 2 || any(r1[diffs] == "*", r2[diffs] == "*"))
    return(NULL)
  
  # Summaries
  summaryAM = summaryPM = NULL
  
  # Case 1: Single vic matching 2 missing
  # TODO: Doesn't handle vic matching 3 or more missing
  if(sum(diffs) == 1) {
    vic = vics[diffs]
    miss = .myintersect(missing, j[1:2, vic]) # intersect: for sorting
    endIdx = match(FALSE, j[[vic]] %in% miss) # 1st row not matching top 2
    GLR = exp(j$loglik[1] - j$loglik[endIdx])
    if(is.na(GLR) || GLR < threshold)
      return(NULL)

    rel = verbalisr::verbalise(dvi$am, miss) |> 
      format(cap = FALSE, simplify = TRUE, collapse = " + ")
    summaryAM = data.frame(Family = fam, Missing = miss, Sample = vic, GLR = GLR, 
                        Conclusion = "Symmetric match", 
                        Comment = sprintf("%s also matches %s (%s)", vic, rev(miss), rel),  
                        row.names = NULL)
    
    summaryPM = data.frame(Sample = vic, Family = fam, Missing = paste(miss, collapse = "/"),
                        GLR = GLR, Conclusion = "Symmetric match", 
                        Comment = sprintf("%s and %s are %s", miss[1], miss[2], rel), 
                        row.names = NULL)
  }
  else if(sum(diffs) == 2 && setequal(r1[diffs], r2[diffs])) {  # Symmetric pair
    vic = vics[diffs]
    miss = .myintersect(missing, r1[diffs]) # intersect: for sorting
    endIdx = match(FALSE, j[[vic[1]]] %in% miss & j[[vic[2]]] %in% miss) # 1st row not matching
    GLR = exp(j$loglik[1] - j$loglik[endIdx])
    if(is.na(GLR) || GLR < threshold)
      return(NULL)

    rel = verbalisr::verbalise(dvi$am, miss) |> 
      format(cap = TRUE, simplify = TRUE, collapse = " + ")
    
    summaryAM = data.frame(Family = fam, Missing = miss, 
                        Sample = paste(vic, collapse = "/"), GLR = GLR, 
                        Conclusion = "Symmetric match", 
                        Comment = sprintf("%s: {%s} = {%s}", rel, toString(miss), toString(vic)),
                        row.names = NULL)
    
    summaryPM = data.frame(Sample = vic, Family = fam, Missing = paste(miss, collapse = "/"),
                        GLR = GLR, Conclusion = "Symmetric match", 
                        Comment = sprintf("%s: {%s} = {%s}", rel, toString(vic), toString(miss)),
                        row.names = NULL)
  }
  
  # Return summaries
  list(AM = summaryAM, PM = summaryPM)
}
