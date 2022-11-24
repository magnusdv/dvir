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
  
  if(!length(missing))
    stop2("The dataset has no missing persons")
  
  if(!length(pm))
    stop2("Empty PM data")
  
  if(!length(am))
    stop2("Empty AM data")
  
  # Enforce lists
  if(is.singleton(pm)) 
    pm = list(pm)
   
  # Ensure `pm` is named
  names(pm) = unlist(labels(pm))
  
  if(!(is.ped(am) || is.pedList(am)))
    stop2("Argument `am` must be a `ped` object or a list of such")
  
  if(is.null(missing))
    stop2("Argument `missing` cannot be NULL")
  
  if(!all(missing %in% unlist(labels(am))))
    stop2("Missing person not found in AM data: ", setdiff(missing, unlist(labels(am))))
  
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
#' @param pm A vector with names or indices of victim samples. By
#'   default, all are included.
#' @param am A vector with names or indices of AM components. By
#'   default, components without remaining missing individuals are dropped.
#' @param missing A vector with names or indices of missing persons. By
#'   default, all missing persons in the remaining AM families are included. 
#' @param verbose A logical.
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
subsetDVI = function(dvi, pm = NULL, am = NULL, missing = NULL, verbose = TRUE) {
  dvi = consolidateDVI(dvi)
  
  pmNew = dvi$pm
  amNew = dvi$am
  missNew = dvi$missing
  
  if(!is.null(pm)) {
    pmNew = dvi$pm[pm]
    
    err = vapply(pmNew, is.null, FALSE)
    if(any(err))
      stop2("Unknown name/index of PM singleton: ", pm[err])
  }
  
  if(!is.null(am)) {
    amNew = dvi$am[am]
    
    err = vapply(amNew, is.null, FALSE)
    if(any(err))
      stop2("Unknown name/index of AM family: ", am[err])
    
    if(is.null(missing)) {
      comps = getComponent(amNew, dvi$missing, checkUnique = FALSE, errorIfUnknown = FALSE)
      if(anyNA(comps)) {
        message("Also removing missing persons in the dropped AM families: ", toString(dvi$missing[!is.na(comps)]))
        missNew = dvi$missing[!is.na(comps)]
      }
    }
  }
  
  if(!is.null(missing)) {
    
    # Add names to allow subsetting by name or index
    miss0 = dvi$missing
    names(miss0) = miss0
  
    missNew = miss0[missing]
    
    err = is.na(missNew)
    if(any(err))
      stop2("Unknown name/index of missing person: ", missing[err])
    
    # If AM subset not given by used, remove components without missing persons
    if(is.null(am)) {
      comps = getComponent(dvi$am, missNew, checkUnique = FALSE, errorIfUnknown = FALSE)
      if(anyNA(comps))
        stop2("Missing person not found in AM data: ", missNew[is.na(comps)])
      
      unused = setdiff(seq_along(dvi$am), comps)
      if(length(unused)) {
        message("Also removing AM families with no missing persons: ", toString(unused))
        amNew[unused] = NULL
      }
    }
  }
  
  dviData(pmNew, amNew, missNew)
}
