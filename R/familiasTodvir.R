#' Convert a Familias file to DVI data
#'
#' This is a wrapper for [readFam()] that reads Familias files with DVI
#' information.
#'
#' @inheritParams relabelDVI
#' @param famfile Path to Familias file.
#' @param verbose A logical. Passed on to [readFam()].
#' @param character missingIdentifier.
#'  First part of names used for the missing
#'
#' @return A `dviData` object.
#'
#' @details The sex of the missing persons need to be checked as this
#'   information may not be correctly recorded in the fam file.
#'
#' @seealso [jointDVI()], [dviData()], [relabelDVI()]
#'
#' @examples
#'
#' # Family with three missing
#' file = "https://familias.name/dviapp/example8.fam"
#' 
#' # Read file without relabelling
#' y = familiasTodvir(file)
#' plotDVI(y)
#' 
#' # With relabelling
#' z = familiasTodvir(file, missingFormat = "M[FAM]-[IDX]",
#'                    refPrefix = "ref", othersPrefix = "E")
#' plotDVI(z)
#'
#' \dontrun{
#' # No data for missing, an error is returned:
#' familiasTodvir("https://familias.name/dviapp/BrotherPower.fam")
#' }
#'
#' @importFrom forrel readFam
#' @export
familiasTodvir = function(famfile, victimPrefix = NULL, familyPrefix = NULL,
                          refPrefix = NULL, missingPrefix = NULL, 
                          missingFormat = NULL, othersPrefix = NULL,
                          verbose = FALSE, missingIdentifier = "^Missing"){
  
  # Return an error if there is no DVI input. Code copied from `readFam`.
  raw = readLines(famfile)
  x = gsub("\\\"", "", raw)
  if (!"[DVI]" %in% x)
    stop2("No DVI input found. Use `forrel::readFam()` to read the file.")
  
  x = readFam(famfile, verbose = verbose)

  if(length(x) == 1) # There are no reference families, no information on the missing.
    stop2("No reference families found. Use `forrel::readFam()` to read the file.")
  

  # PM data -----------------------------------------------------------------
  
  pm = x$`Unidentified persons`
  

  # AM data -----------------------------------------------------------------
  am = lapply(x[-1], function(dat) {
      ref = dat[[2]]
      if(!is.ped(ref))
        ref = ref[[which(pedsize(ref) > 1)]]     
      ref
    })

  
  # Check for identically named reference individuals
  if(!is.null(am)){
    if(any(duplicated(typedMembers(am))))
      stop2("Typed members of reference families must be named differently.")
  }
  
  
  # Missing individuals -----------------------------------------------------

  missing = grep(missingIdentifier, unlist(labels(am)), value = TRUE)
  
  if(!length(missing))
    stop2("Reference pedigree without missing")
  

  # Create DVI object -------------------------------------------------------

  dvi0 = dviData(pm = pm, am = am, missing = missing)
  
  # Relabel missing persons
  dvi = relabelDVI(dvi0, victimPrefix = victimPrefix, familyPrefix = familyPrefix,
                   refPrefix = refPrefix, missingPrefix = missingPrefix, 
                   missingFormat = missingFormat, othersPrefix = othersPrefix)
  dvi
}