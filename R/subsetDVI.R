#' Extract a subset of a DVI dataset
#'
#' @param dvi A [dviData()] object
#' @param pm A vector with names or indices of victim samples. By default, all
#'   are included.
#' @param am A vector with names or indices of AM components. By default,
#'   components without remaining missing individuals are dropped.
#' @param missing A vector with names or indices of missing persons. By default,
#'   all missing persons in the remaining AM families are included.
#' @param removeUnpairedPM A logical, by default TRUE, removing PM samples with
#'   no remaining pairings.
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
subsetDVI = function(dvi, pm = NULL, am = NULL, missing = NULL, removeUnpairedPM = TRUE, verbose = TRUE) {
  
  if(verbose)
    cat("Reducing DVI dataset\n")
  
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
      
      if(!length(amNew))
        missNew = character(0)
      else {
        comp = getComponent(amNew, dvi$missing, checkUnique = FALSE, errorIfUnknown = FALSE)
        missNew = dvi$missing[!is.na(comp)]
      }
      missRem = .mysetdiff(dvi$missing, missNew)
      if(verbose && length(missRem))
        cat(sprintf("Removing %s in excluded families: %s\n", 
                    pluralise("missing person", length(missRem)), trunc(missRem, 6)))
    }
  }
  
  if(!is.null(missing)) {
    
    # Add names to allow subsetting by name or index
    miss0 = dvi$missing
    names(miss0) = miss0
    
    missNew = unname(miss0[missing])
    
    err = is.na(missNew)
    if(any(err))
      stop2("Unknown name/index of missing person: ", missing[err])
    
    # If AM subset not given by user, remove components without missing persons
    if(is.null(am)) {
      comps = getComponent(dvi$am, missNew, checkUnique = FALSE, errorIfUnknown = FALSE)
      if(anyNA(comps))
        stop2("Missing person not found in AM data: ", missNew[is.na(comps)])
      
      unused = setdiff(seq_along(dvi$am), comps)
      if(length(unused)) {
        if(verbose) {
          nun = length(unused)
          famNames = names(amNew)[unused] %||% unused
          cat(sprintf("Removing %d AM famil%s with no remaining missing persons: %s\n", 
                          nun, if(nun == 1) "y" else "ies", trunc(famNames, 6)))
        }
        amNew[unused] = NULL
      }
    }
  }
  
  dviNew = dviData(pmNew, amNew, missNew)
  
  # Fix pairings if they were included in original dataset
  if(!is.null(dvi$pairings)) {
    removedMiss = setdiff(dvi$missing, missNew)
    dviNew$pairings = lapply(dvi$pairings[names(pmNew)], function(v) v[!v %in% removedMiss])
  } 
    
  # PMs with no remaining pairings?
  if(removeUnpairedPM) {
    excl = vapply(dviNew$pairings, function(v) length(v) == 1 && v == "*", FUN.VALUE = FALSE)
    if(sum(excl) > 0) {
      if(verbose)
        cat(sprintf("Removing %s with no remaining pairings: %s\n", 
                    pluralise("PM sample", sum(excl)), trunc(names(excl)[excl], 6)))
      dviNew$pm[excl] = dviNew$pairings[excl] = NULL
    }
  }
  
  dviNew
}
