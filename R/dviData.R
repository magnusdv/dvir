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

