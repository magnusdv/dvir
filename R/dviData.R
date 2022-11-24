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
  
  # Enforce lists
  if(is.singleton(pm)) 
    pm = list(pm)
   
  # Ensure `pm` is named
  names(pm) = unlist(labels(pm))
  
  if(is.null(am) || !(is.ped(am) || is.pedList(am)))
    stop2("Argument `am` must be a `ped` object or a list of such")
  
  if(is.null(missing))
    stop2("Argument `missing` cannot be NULL")
  
  if(anyDuplicated(missing))
    stop2("Duplicated entries of `missing`: ", missing[duplicated(missing)])
  
  amlist = if(is.ped(am)) list(am) else am
  comps = getComponent(amlist, missing, checkUnique = TRUE, errorIfUnknown = FALSE)
  if(anyNA(comps))
    stop2("Missing person not found in AM data: ", missing[is.na(comps)])
  
  unused = setdiff(seq_along(amlist), comps)
  if(length(unused))
    warning("Some components of `am` have no missing individuals: ", unused)
  
  structure(list(pm = pm, am = am, missing = missing), 
                  class = "dviData")
}


#' @export
print.dviData = function(x, ..., heading = "DVI dataset:", printMax = 10) {
  dvi = consolidateDVI(x)
  pm = dvi$pm
  am = dvi$am
  missing = dvi$missing
  
  vics = names(pm)
  refs = typedMembers(am)
  nam = length(am)
  
  message(heading)
  message(sprintf(" %d victims: %s", length(pm), trunc(vics, printMax)))
  message(sprintf(" %d missing: %s", length(missing), trunc(missing, printMax)))
  message(sprintf(" %d typed refs: %s", length(refs), trunc(refs, printMax)))
  message(sprintf(" %d reference famil%s", nam, ifelse(nam == 1, "y", "ies")))
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



#' Extract a subset of a DVI dataset
#'
#' @param dvi A [dviData()] object
#' @param pm A vector with names or indices of victim samples to include. By
#'   default, all are included.
#' @param am A vector with names or indices of AM components include. By
#'   default, all relevant components are included.
#' @param missing A vector with names or indices of missing persons to include. By
#'   default, all relevant missing persons are included. 
#'
#' @return A `dviData` object.
#'
#' @examples
#' 
#' subsetDVI(example2, pm = 1:2) |> plotDVI()
#' subsetDVI(example2, pm = "V1", am = 1) |> plotDVI()
#' subsetDVI(example2, missing = "M3") |> plotDVI()
#' 
#' @export
subsetDVI = function(dvi, pm = NULL, am = NULL, missing = NULL) {
  dvi = consolidateDVI(dvi)
  
  pmNew = dvi$pm
  amNew = dvi$am
  missNew = dvi$missing
  
  if(!is.null(pm)) 
    pmNew = dvi$pm[pm]
  
  if(!is.null(am)) {
    amNew = dvi$am[am]
    
    if(is.null(missing)) {
      comps = getComponent(amNew, dvi$missing, checkUnique = FALSE, errorIfUnknown = FALSE)
      missNew = dvi$missing[!is.na(comps)]
    }
  }
  
  if(!is.null(missing)) {
    
    if(is.character(missing))
      missNew = missing
    else
      missNew = dvi$missing[missing]
    
    # If AM subset not given by used, remove components without missing persons
    if(is.null(am)) {
      comps = getComponent(dvi$am, missNew, checkUnique = FALSE, errorIfUnknown = FALSE)
      if(anyNA(comps))
        stop2("Missing person not found in AM data: ", missNew[is.na(comps)])
      
      amNew = dvi$am[unique.default(comps)]
    }
  }
  
  dviData(pmNew, amNew, missNew)
}
