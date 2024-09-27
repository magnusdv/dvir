#' DVI data
#'
#' @param pm A list of singletons: The victim samples.
#' @param am A list of pedigrees: The reference families.
#' @param missing A character vector with names of missing persons.
#' @param generatePairings A logical. If TRUE (default) a list of sex-compatible
#'   pairings is included as part of the output.
#' @param ignoreSex A logical.
#' @param dvi A `dviData` object.
#' @param pairings A list of pairings.
#' @param errorIfEmpty A logical.
#' @param verbose A logical.
#'
#' @return An object of class `dviData`, which is basically a list of `pm`,
#'   `am`, `missing` and `pairings`.
#'
#' @examples
#' dvi = dviData(pm = singleton("V1"), am = nuclearPed(1), missing = "3")
#' dvi
#' 
#' checkDVI(dvi)
#' 
#' @export
dviData = function(pm, am, missing, generatePairings = TRUE) {
  if(inherits(pm, "dviData"))
    return(pm)
  
  # Enforce lists
  if(is.singleton(pm)) 
    pm = list(pm)
   
  # Ensure `pm` is named
  names(pm) = unlist(labels(pm))
  
  if(!(length(am) == 0 || is.ped(am) || is.pedList(am)))
    stop2("Argument `am` must be a `ped` object or a list of such")
  
  if(is.null(missing))
    stop2("Argument `missing` cannot be NULL")
  
  if(!all(missing %in% unlist(labels(am))))
    stop2("Missing person not found in AM data: ", setdiff(missing, unlist(labels(am))))
  
  missing = as.character(missing)
  
  dvi = structure(list(pm = pm, am = am, missing = missing, pairings = NULL), 
                  class = "dviData")
  
  if(generatePairings) 
    dvi$pairings = generatePairings(dvi, ignoreSex = FALSE)
  
  consolidateDVI(dvi)
}


#' @export
print.dviData = function(x, ..., heading = "DVI dataset:", printMax = 10) {
  
  dvi = x
  pm = dvi$pm
  am = dvi$am
  missing = dvi$missing
  
  npm = length(pm)
  nam = length(am)
  
  vics = names(pm)
  refs = if(nam) typedMembers(am) else NULL
  amNames = names(am)
  
  # Number of victim males and females:
  nVfemales = length(females(pm))
  nVmales = length(males(pm))
  
  # Number of missing males and females:
  nMPsex = tabulate(getSex(x$am, missing), nbins = 2)
  
  # Need to know later if there are no, equal (>0) or different number of markers
  rangePM = if(npm) range(nMarkers(pm, compwise = TRUE)) else c(0,0)
  rangeAM = if(nam) range(nMarkers(am, compwise = TRUE)) else c(0,0)
  
    
  cat(heading, "\n")
  cat(sprintf(" %d victim%s (%dM/%dF): %s\n", 
              length(pm), if(length(pm) == 1) "" else "s", nVmales, nVfemales, trunc(vics, printMax)))
  cat(sprintf(" %d missing (%dM/%dF): %s\n", 
              length(missing), nMPsex[1], nMPsex[2], trunc(missing, printMax)))
  cat(sprintf(" %d typed ref%s: %s\n", length(refs), if(length(refs) == 1) "" else "s", trunc(refs, printMax)))
  cat(sprintf(" %d ref famil%s: %s\n", 
              nam, ifelse(nam == 1, "y", "ies"), trunc(amNames, printMax)))
  
  ### Number of markers
  
  # Simple case: PM and AM equal, same number for all
  if(min(rangePM, rangeAM) == max(rangePM, rangeAM))
    cat("Number of markers, PM and AM:", rangePM[1], "\n")
  else {
    if(npm > 0) {
      if(rangePM[1] == rangePM[2])
        cat("Number of markers, PM:", rangePM[1], "\n")
      else 
        cat(sprintf("Number of markers, PM: Ranges from %d to %d\n", 
                    rangePM[1], rangePM[2]))
    }
    if(nam > 0) {
      if(rangeAM[1] == rangeAM[2])
        cat("Number of markers, AM:", rangeAM[1], "\n")
      else 
        cat(sprintf("Number of markers, AM: Ranges from %d to %d\n", 
                    rangeAM[1], rangeAM[2]))
    }
  }
}


# Utility to be run in the beginning of all major functions, 
# ensuring that the input is correctly formatted.
# Particularly important if the user has modified the `dvi` object.
consolidateDVI = function(dvi, dedup = FALSE) {

  if(!inherits(dvi, "dviData"))
    stop2("Cannot consolidate; input is not a `dviData` object")
  
  # Make sure pm is a list
  if(is.ped(dvi$pm))   # include non-singletons, so that they can be caught below
    dvi$pm = list(dvi$pm)
  
  # Ensure `pm` is correctly named
  labs = unlist(lapply(dvi$pm, function(x) x$ID), use.names = FALSE) # faster than labels(pm)
  if(length(labs) != length(dvi$pm))
    stop2("PM part is not a list of singletons")
  
  names(dvi$pm) = labs
  
  # Make sure am is a list
  if(is.ped(dvi$am))
    dvi$am = list(dvi$am)
  
  # Enforce family names if not present (F1, F2, ...)
  if(is.null(names(dvi$am)) && length(dvi$am))
    names(dvi$am) = paste0("F", seq_along(dvi$am))

  # Ensure `missing` is an unnamed character
  dvi$missing = as.character(dvi$missing)
  
  if(dedup)
    dvi = relabelDVI(dvi, othersPrefix = "")
  
  dvi
}



#' @rdname dviData
#' @export
checkDVI = function(dvi, pairings = NULL, errorIfEmpty = FALSE, 
                    ignoreSex = FALSE, verbose = TRUE){
  
  #TODO: Check for duplicated labels! (e.g. among refs)
  
  if(verbose)
    cat("Checking DVI dataset consistency\n")
  
  # Assumes `dvi` has been consolidated
  
  pm = dvi$pm
  am = dvi$am
  missing = dvi$missing
  
  # Added to avoid crash in certain cases.
  if(length(pm) == 0 || length(missing) == 0) {
    if(errorIfEmpty) stop2("Empty DVI problem") 
    else return()
  }
  
  # Check that PM is a list of singletons
  if(!is.list(pm) || !all(vapply(pm, is.singleton, TRUE)))
    stop2("PM data is not a list of singletons")
  
  # Check that all missing are members of a ref pedigree
  comps = getComponent(am, missing, checkUnique = TRUE, errorIfUnknown = FALSE)
  if(anyNA(comps))
    stop2("Missing person not found in the AM pedigree(s): ", missing[is.na(comps)])
  
  unused = setdiff(seq_along(am), comps)
  if(length(unused))
    warning("AM families with no missing individuals: ", toString(unused), 
            call. = FALSE, immediate. = TRUE)
  
  # Check if marker sets are identical
  markersAM = name(am)
  markersPM = name(pm)
  if(!setequal(markersAM, markersPM)) {
    msg = "Warning: Unequal marker sets in PM and AM\n"
    if(length(notinPM <- setdiff(markersAM, markersPM))) {
      n1 = length(notinPM)
      line1 = sprintf(" %d marker%s in AM but not in PM: %s\n", n1, if(n1 == 1) "" else "s", toString(notinPM))
      msg = paste0(msg, line1)
    }
    if(length(notinAM <- setdiff(markersPM, markersAM))) {
      n2 = length(notinAM)
      line2 = sprintf(" %d marker%s in PM but not in AM: %s\n", n2, if(n2 == 1) "" else "s", toString(notinAM))
      msg = paste0(msg, line2)
    }
    cat(msg)
  }
  
  pairings = pairings %||% dvi$pairings
  
  if(length(pairings) == 0)
    return(invisible()) #stop2("Argument `pairings` has length 0")
  
  vics = names(pm)
  vicSex = getSex(pm, vics, named = TRUE)
  
  candidMP = setdiff(unlist(pairings), "*")
  candidSex = getSex(am, candidMP, named = TRUE)
  
  if(!all(candidMP %in% missing))
    stop2("Pairing error: Candidate is not a missing person: ", setdiff(candidMP, missing))
  
  for(v in vics) {
    candid = pairings[[v]]
    if(length(candid) == 0)
      stop2("No available candidate for victim ", v)
    
    if(any(duplicated(candid)))
      stop2("Duplicated candidate for victim ", v)
    
    cand = setdiff(candid, "*")
    if(length(cand) == 0)
      next
    
    if(!ignoreSex) {
      correctSex = candidSex[cand] == vicSex[v]
      if(!all(correctSex)) 
        stop2("Candidate for victim ", v, " has wrong sex: ", cand[correctSex])
    }
  }
  
  if(verbose)
    cat("No problems found\n")
}



#' Get AM component of selected individuals
#'
#' @param dvi A [dviData()] object.
#' @param ids A vector of ID labels of members of `dvi$am`.
#'
#' @return A vector of the same length as `ids`, containing the family names of
#'   the `ids` individuals.
#'
#' @examples
#' getFamily(example2, ids = example2$missing)
#'
#' @export
getFamily = function(dvi, ids) {
  if(!length(ids))
    return(character(0))
  dvi = consolidateDVI(dvi)
  famnames = names(dvi$am)
  if(is.null(famnames))
    stop2("AM data has no family names")
  comp = getComponent(dvi$am, ids, checkUnique = TRUE, errorIfUnknown = TRUE)
  comp = famnames[comp]
  names(comp) = ids
  comp
}


#' Find the simple families of a DVI dataset
#'
#' Extract the names (if present) or indices of the *simple* reference families,
#' i.e., the families containing exactly 1 missing person.
#'
#' @param dvi A `dviData` object.
#'
#' @return A character (if `dvi$am` has names) or integer vector.
#' @seealso [getFamily()]
#' 
#' @examples
#' # No simple families
#' simple1 = getSimpleFams(example1)
#' stopifnot(length(simple1) == 0)
#' 
#' # Second family is simple
#' simple2 = getSimpleFams(example2)
#' stopifnot(simple2 == "F2")
#' 
#' @export
getSimpleFams = function(dvi) {
  dvi = consolidateDVI(dvi)
  famnames = names(dvi$am)
  
  # AM component of each missing
  fams = getFamily(dvi, ids = dvi$missing)
  
  # Number of missing in each (better than `table`)
  nMiss = tabulate(match(fams, famnames))

  res = famnames[nMiss == 1]
  
  # Convert to integer if indices
  if(is.integer(fams)) 
    res = as.integer(res)
  
  res
}

dviEqual = function(dvi1, dvi2) {
  test1 = identical(dvi1$pm, dvi2$pm) && 
    identical(dvi1$am, dvi2$am) && 
    identical(dvi1$missing, dvi2$missing) &&
    identical(lengths(dvi1$pairings), lengths(dvi2$pairings))
  if(!test1)
    return(FALSE)
  
  # Pairings may come in different order
  
  for(v in names(dvi1$pairings))
    if(!setequal(dvi1$pairings[[v]], dvi2$pairings[[v]]))
      return(FALSE)
  
  return(TRUE)
}

# How many pairings are removed in the reduced dataset?
# Only count those between vics and missing persons still included
removedPairings = function(dviRed, dvi) {
  newvics = names(dviRed$pm)
  
  # Previous pairings for the current vics
  pOld = dvi$pairings[newvics]
  
  # Remove pairings to missing individuals no longer included
  remMiss = .mysetdiff(dvi$missing, dviRed$missing)
  if(length(remMiss))
    pOld = lapply(pOld, function(a) .mysetdiff(a, remMiss))
  
  sum(lengths(pOld)) - sum(lengths(dviRed$pairings))
}
