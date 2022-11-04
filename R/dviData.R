#' DVI data
#'
#'
#' @param pm A list of singletons: The victim samples.
#' @param am A list of pedigrees: The reference families.
#' @param missing A character vector with names of missing persons.
#' @param validate A logical, by default TRUE.
#'
#' @return An object of class `dviData`, which is basically a list of `pm`, `am`
#'   and `missing` after consolidation.
#'
#' @examples
#' dviData(pm = singleton("V1"), am = nuclearPed(1), missing = "3")

#' @export
dviData = function(pm, am, missing, validate = TRUE) {
  if(inherits(pm, "dviData"))
    return(pm)
  
  # Enforce lists
  if(is.singleton(pm)) 
    pm = list(pm)
  if(is.ped(am)) 
    am = list(am)
  
  # Ensure `pm` is named
  names(pm) = unlist(labels(pm))
  
  dvi = structure(list(pm = pm, am = am, missing = missing), 
                  class = "dviData")
  # if(validate)
  #   dvi = validateDvi(dvi)
  
  dvi
}

#' @export
print.dviData = function(x, ..., method = NULL, printMax = 10)
  summariseDVI(x, method = method, printMax = printMax)


#' Plot a DVI problem
#'
#' @param dvi A `dviData` object, typically created with [dviData()].
#' @param hatched A character vector of ID labels, or the name of a function. By
#'   default, typed individuals are hatched.
#' @param col A list of colour vectors (see [pedtools::plot.ped()]). By default,
#'   missing members of `dvi$am` are shown in red.
#' @param ... Further parameters to be passed on to [pedtools::plotPedList()].
#'   Useful parameters include `frames`, `marker`, `newdev` and `titles`.
#'
#' @return NULL
#'
#' @examples
#'
#' plotDVI(example1)
#' 
#' plotDVI(example2, new = T, frames = F, marker = 1, cex = 1.2)
#' 
#' plotDVI(icmp)
#' 
#' @export
plotDVI = function(dvi, hatched = typedMembers, col = list(red = dvi$missing), ...) {
  plotPedList(list(PM = dvi$pm, AM = dvi$am), hatched = hatched, col = col, ...)
}