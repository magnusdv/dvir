#' Convert a Familias file to DVI data
#'
#' This is a wrapper for [pedFamilias::readFam()] that reads Familias files with
#' DVI information.
#'
#' @inheritParams relabelDVI
#' @param famfile Path to Familias file.
#' @param verbose A logical, passed on to `readFam()`.
#' @param missingIdentifier A character of length 1 used to identify missing
#'   persons in the Familias file. The default chooses everyone whose label
#'   begins with "Missing".
#'
#' @return A `dviData` object.
#'
#' @details The sex of the missing persons need to be checked as this
#'   information may not be correctly recorded in the fam file.
#'
#' @seealso [dviData()], [relabelDVI()]
#'
#' @examples
#'
#' # Family with three missing
#' file = system.file("extdata", "dvi-example.fam", package="dvir")
#'
#' # Read file without relabelling
#' y = familias2dvir(file)
#' plotDVI(y)
#'
#' # With relabelling
#' z = familias2dvir(file, missingFormat = "M[FAM]-[IDX]",
#'                    refPrefix = "ref", othersPrefix = "E")
#' plotDVI(z)
#'
#' @export
familias2dvir = function(famfile, victimPrefix = NULL, familyPrefix = NULL,
                         refPrefix = NULL, missingPrefix = NULL, 
                         missingFormat = NULL, othersPrefix = NULL,
                         verbose = FALSE, missingIdentifier = "^Missing"){
  
  # Read fam file
  x = pedFamilias::readFam(famfile, useDVI = TRUE, verbose = verbose)

  # PM data -----------------------------------------------------------------
  
  pm = x$`Unidentified persons`
  if(is.null(pm))
    stop2("No `Unidentified persons` found")

  # AM data -----------------------------------------------------------------
  
  am = x[-1]
  
  # Remove untyped components
  am = lapply(names(am), function(refnm) {
    ref = am[[refnm]]
    if(is.ped(ref))
      return(ref)
    
    # Remove Reference pedigree if present
    idx = match("Reference pedigree", names(ref), nomatch = 0)
    if(length(ref) == 2 && idx > 0)
      ref = ref[[3 - idx]] # the other
    
    # Expect single component with typed references
    cmp = getComponent(ref, typedMembers(ref))
    if(max(cmp) > min(cmp))
      stop2("Disconnected reference family: ", refnm)
    
    ref[[cmp[1]]]
  })

  # Check for duplicated names among reference individuals
  if(!is.null(am)) {
    typed = typedMembers(am)
    dups = anyDuplicated.default(typed)
    if(dups)
      stop2("Duplicated name among reference individuals: ", typed[dups])
  }
  
  # Missing individuals -----------------------------------------------------

  missing = grep(missingIdentifier, unlist(labels(am)), value = TRUE)
  
  if(!length(missing))
    stop2("Reference pedigree without missing")
  

  # Create DVI object -------------------------------------------------------

  dvi0 = dviData(pm = pm, am = am, missing = missing, generatePairings = FALSE)
  
  # Relabel missing persons (NB: pairings are generated here)
  relabelDVI(dvi0, victimPrefix = victimPrefix, familyPrefix = familyPrefix,
             refPrefix = refPrefix, missingPrefix = missingPrefix, 
             missingFormat = missingFormat, othersPrefix = othersPrefix)
}