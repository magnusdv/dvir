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



consolidate = function(dvi) {
  if(is.ped(dvi$am))
    dvi$am = list(dvi$am)
  
  dvi
}
