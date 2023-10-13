#' DVI data
#'
#' @param pm A list of singletons: The victim samples.
#' @param am A list of pedigrees: The reference families.
#' @param missing A character vector with names of missing persons.
#' @param generatePairings A logical. If TRUE (default) a list of sex-compatible
#'   pairings is included as part of the output.
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
  
  dvi = structure(list(pm = pm, am = am, missing = missing, pairings = NULL), 
                  class = "dviData")
  
  if(generatePairings) 
    dvi$pairings = generatePairings(dvi, ignoreSex = FALSE)
  
  consolidateDVI(dvi)
}


#' @export
print.dviData = function(x, ..., heading = "DVI dataset:", printMax = 10) {
  dvi = consolidateDVI(x)
  pm = dvi$pm
  am = dvi$am
  missing = dvi$missing
  
  npm = length(pm)
  nam = length(am)
  
  vics = names(pm)
  refs = if(nam) typedMembers(am) else NULL
  amNames = names(am) %||% "(unnamed)"
  
  # Number of victim males and females:
  nVfemales = length(females(pm))
  nVmales = length(males(pm))
  
  # Number of missing males and females:
  nMPsex = tabulate(getSex(x$am, missing), nbins = 2)
  
  # Need to know later if there are no, equal (>0) or different number of markers
  nMarkersPM = if(npm) range(nMarkers(pm, compwise = TRUE)) else 0
  nMarkersAM = if(nam) range(nMarkers(am, compwise = TRUE)) else 0
  
    
  cat(heading, "\n")
  cat(sprintf(" %d victims (%dM/%dF): %s\n", 
              length(pm), nVmales, nVfemales, trunc(vics, printMax)))
  cat(sprintf(" %d missing (%dM/%dF): %s\n", 
              length(missing), nMPsex[1], nMPsex[2], trunc(missing, printMax)))
  cat(sprintf(" %d typed refs: %s\n", length(refs), trunc(refs, printMax)))
  cat(sprintf(" %d ref famil%s: %s\n", 
              nam, ifelse(nam == 1, "y", "ies"), trunc(amNames, printMax)))
  
  ### Number of markers
  
  # Simple case: PM and AM equal, same number for all
  if(min(nMarkersPM, nMarkersAM) == max(nMarkersPM, nMarkersAM))
    cat("Number of markers, PM and AM:", nMarkersPM[1], "\n")
  else {
    if(nMarkersPM[1] == nMarkersPM[2])
      cat("Number of markers, PM:", nMarkersPM[1], "\n")
    else 
      cat(sprintf("Number of markers, PM: Ranges from %d to %d\n", 
                  nMarkersPM[1], nMarkersPM[2]))
    if(nMarkersAM[1] == nMarkersAM[2])
      cat("Number of markers, AM:", nMarkersAM[1], "\n")
    else 
      cat(sprintf("Number of markers, AM: Ranges from %d to %d\n", 
                  nMarkersAM[1], nMarkersAM[2]))
  }
}


# Utility to be run in the beginning of all major functions, 
# ensuring that the input is correctly formatted.
# Particularly important if the user has modified the `dvi` object.
consolidateDVI = function(dvi) {
  
  if(!inherits(dvi, "dviData"))
    stop2("Cannot consolidate; input is not a `dviData` object")
  
  # Make sure pm is a list
  if(is.singleton(dvi$pm)) 
    dvi$pm = list(dvi$pm)
  
  # Ensure `pm` is correctly named
  labs = unlist(lapply(dvi$pm, function(x) x$ID), use.names = FALSE) # faster than labels(pm)
  if(length(labs) != length(dvi$pm))
    stop2("PM part is not a list of singletons")
  
  names(dvi$pm) = labs
  
  # Make sure am is a list
  if(is.ped(dvi$am))
    dvi$am = list(dvi$am)
  
  dvi
}



#' @rdname dviData
#' @export
checkDVI = function(dvi, pairings = NULL, errorIfEmpty = FALSE, ignoreSex = FALSE){
  
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
    msg = "Unequal marker sets in PM and AM:"
    if(length(notinPM <- setdiff(markersAM, markersPM)))
      msg = paste0(msg, "\n * marker(s) in AM but not in PM: ", toString(notinPM))
    if(length(notinAM <- setdiff(markersPM, markersAM)))
      msg = paste0(msg, "\n * marker(s) in PM but not in AM: ", toString(notinAM))
    warning(msg, call. = FALSE, immediate. = TRUE)
  }
  
  if(is.null(pairings))
    return(invisible())
  
  if(length(pairings) == 0)
    stop2("Argument `pairings` has length 0")
  
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
}



#' Get AM component of selected individuals
#'
#' @param dvi A [dviData()] object.
#' @param ids A vector of ID labels of members of `dvi$am`.
#'
#' @return A vector of the same length as `ids`, containing the family names (if
#'   `dvi$am` is named) or component indices (otherwise) of the `ids`
#'   individuals.
#'
#' @examples
#' getFamily(example2, ids = example2$missing)
#' 
#' @export
getFamily = function(dvi, ids) {
  comp = getComponent(dvi$am, ids, checkUnique = TRUE, errorIfUnknown = TRUE)
  if(!is.null(famnames <- names(dvi$am)))
    comp = famnames[comp]
  names(comp) = ids
  comp
}