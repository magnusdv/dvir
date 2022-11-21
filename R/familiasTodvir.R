#' Convert familias file to DVI data 
#'
#' This is a wrapper for [readFam()] that reads Familias files with DVI
#' information.
#'
#' @param famfile Path to Familias file.
#' @param victimPrefix Prefix used to label PM individuals.
#' @param familyPrefix Prefix used to label the AM families.
#' @param refPrefix Prefix used to label the reference individuals, i.e., the
#'   typed members of the AM families.
#' @param missingPrefix Prefix used to label the missing persons in the AM
#'   families. The word "family" is treated as a special case, where the family
#'   name is used as prefix in each family, e.g., F1-1, F1-2, F2-1, ...
#' @param missingFormat A string indicating family-wise labelling of missing
#'   persons, using `[FAM]` an `[IDX]` as place holders for the family index and
#'   the missing person index within the family. See Examples in [relabelDVI()].
#' @param othersPrefix	Prefix used to label other untyped individuals. Default: 1, 2, ...

#' @param verbose A logical. Passed on to [readFam()].
#'
#' @return A `dviData` object. 
#'
#' @details The sex of the missing persons need to be checked as this
#' information may not be correctly recorded in the fam file.  
#' At most one of `missingPrefix` and `missingFormat` can be used.
#'
#' @seealso [jointDVI()], [dviData()], [relabelDVI()]
#' @export
#'
#' @examples
#' 
#' # Family with three missing
#' file = "https://familias.name/dviapp/example8.fam"
#' familiasTodvir(file)
#' z = familiasTodvir(file, victimPrefix = "vic", familyPrefix = "fam", 
#'                    refPrefix = "ref", missingPrefix = NULL, 
#'                    missingFormat = "M[FAM]-[IDX]", othersPrefix = "E")
#' z$am[1]
#'
#' \dontrun{
#' # No data for missing, an error is returned:
#' familiasTodvir("https://familias.name/dviapp/BrotherPower.fam")
#' }
#'
#' @importFrom forrel readFam
#' 
familiasTodvir = function(famfile, victimPrefix = NULL, familyPrefix = NULL,
                          refPrefix = NULL, missingPrefix = NULL, 
                          missingFormat = NULL, othersPrefix = NULL,
                          verbose = FALSE){
  
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

  missing = grep("^Missing", unlist(labels(am)), value = TRUE)
  
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