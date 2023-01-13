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
  
  if(verbose)
    message("Reducing DVI dataset")
  
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
        if(verbose)
          message("Removing missing persons in dropped AM families:\n ", 
                  toString(dvi$missing[!is.na(comps)]))
        missNew = dvi$missing[!is.na(comps)]
      }
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
    
    # If AM subset not given by used, remove components without missing persons
    if(is.null(am)) {
      comps = getComponent(dvi$am, missNew, checkUnique = FALSE, errorIfUnknown = FALSE)
      if(anyNA(comps))
        stop2("Missing person not found in AM data: ", missNew[is.na(comps)])
      
      unused = setdiff(seq_along(dvi$am), comps)
      if(length(unused)) {
        if(verbose)
          message("Removing AM families with no missing persons: ", toString(unused))
        amNew[unused] = NULL
      }
    }
  }
  
  dviData(pmNew, amNew, missNew)
}
