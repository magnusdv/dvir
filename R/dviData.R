#' DVI data
#'
#' @param pm A list of singletons: The victim samples.
#' @param am A list of pedigrees: The reference families.
#' @param missing A character vector with names of missing persons.
#'
#' @return An object of class `dviData`, which is basically a list of `pm`, `am`
#'   and `missing`.
#'
#' @examples
#' dviData(pm = singleton("V1"), am = nuclearPed(1), missing = "3")
#' 
#' @export
dviData = function(pm, am, missing) {
  if(inherits(pm, "dviData"))
    return(pm)
  
  #if(!length(missing))
  #  stop2("The dataset has no missing persons")
  
  #if(!length(pm))
  #  stop2("Empty PM data")
  
  #if(!length(am))
  #  stop2("Empty AM data")
  
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
  
  dvi = structure(list(pm = pm, am = am, missing = missing), 
                  class = "dviData")
  
  dvi$pairings = generatePairings(dvi, ignoreSex = FALSE)
  dvi
}


#' @export
print.dviData = function(x, ..., heading = "DVI dataset:", printMax = 10) {
  dvi = consolidateDVI(x)
  pm = dvi$pm
  am = dvi$am
  missing = dvi$missing
  
  vics = names(pm)
  refs = if(length(am)) typedMembers(am) else NULL
  nam = length(am)
  amNames = names(am) %||% "(unnamed)"
  message(heading)
  message(sprintf(" %d victims: %s", length(pm), trunc(vics, printMax)))
  message(sprintf(" %d missing: %s", length(missing), trunc(missing, printMax)))
  message(sprintf(" %d typed refs: %s", length(refs), trunc(refs, printMax)))
  message(sprintf(" %d ref famil%s: %s", 
                  nam, ifelse(nam == 1, "y", "ies"), trunc(amNames, printMax)))
}


# Utility to be run in the beginning of all major functions, 
# ensuring that the input is correctly formatted.
# Particularly important if the user has modified the `dvi` object.
consolidateDVI = function(dvi) {
  
  if(!inherits(dvi, "dviData")) {
    if(setequal(names(dvi), c("pm", "am", "missing")))
      dvi = dviData(pm = dvi$pm, am = dvi$am, missing = dvi$missing)
    else
      stop2("Cannot consolidate `dviData`: Wrong entry names")
  }
  
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



# @rdname jointDVI
# @export
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
    stop2("`pm` object is not a list of singletons")
  
  # Check that all missing are members of a ref pedigree
  comps = getComponent(am, missing, checkUnique = TRUE, errorIfUnknown = FALSE)
  if(anyNA(comps))
    stop2("Missing person not found in the AM pedigree(s): ", missing[is.na(comps)])
  
  unused = setdiff(seq_along(am), comps)
  if(length(unused))
    warning("Some components of `am` have no missing individuals: ", toString(unused), 
            call. = FALSE, immediate. = TRUE)
  
  if(is.null(pairings))
    return(invisible())
  if(length(pairings) == 0)
    stop2("Argument `pairings` has length 0")
  
  vics = names(pm)
  vicSex = getSex(pm, vics, named = TRUE)
  
  candidMP = setdiff(unlist(pairings), "*")
  candidSex = getSex(am, candidMP, named = TRUE)
  
  if(!all(candidMP %in% missing))
    stop2("Indicated pairing candidate is not a missing person: ", setdiff(candidMP, missing))
  
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


#' Summarise a DVI problem
#'
#' DEPRECATED: USE `print(dvi)` INSTEAD. Prints a summary of a given DVI problem, including the number of victims,
#' missing persons, reference families and typed reference individuals. NOTE:
#' 
#' @param dvi A `dviData` object, typically created with [dviData()].
#' @param method A character, used by other methods.
#' @param printMax A positive integer. Vectors longer than this are truncated.
#'
#' @return No return value, called for side effects.
#'
#' @examples
#' summariseDVI(planecrash)
#' summariseDVI(planecrash, printMax = 5)
#'
#' @export
summariseDVI = function(dvi, method = NULL, printMax = 10) {
  
  cat("Note: `summariseDVI()` has been deprecated and merged with the print method for `dviData` objects.\n\n")
  
  # Ensure proper dviData object
  dvi = consolidateDVI(dvi)
  
  pm = dvi$pm
  am = dvi$am
  missing = dvi$missing
  
  vics = names(pm)
  refs = typedMembers(am)
  nam = length(am)
  
  message("DVI problem:")
  message(sprintf(" %d victims: %s", length(pm), trunc(vics, printMax)))
  message(sprintf(" %d missing: %s", length(missing), trunc(missing, printMax)))
  message(sprintf(" %d typed refs: %s", length(refs), trunc(refs, printMax)))
  message(sprintf(" %d reference famil%s", nam, ifelse(nam == 1, "y", "ies")))
  if(!is.null(method))
    message("\n", method)
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