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


# Utility to be run in the beginning of all major functions, ensuring that the input is correctly formatted.
# This is particularly important if the user has modified the `dvi` object.
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
