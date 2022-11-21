#' Convert familias file to dvir input
#'
#' This is a wrapper for [readFam()] that reads Familias files with DVI
#' information.
#'
#' @param famfile Path to Familias file.
#' @param missingPrefix A string.
#' @param verbose A logical. Passed onto [readFam()]
#'
#' @return A `dviData` object. The missing persons are renamed. The first number
#'   indicates the family, the second the missing person in the family.
#'
#' @details The sex of the missing persons need to be checked as this
#'   information may not be correctly recorded in the fam file.
#'
#' @seealso [jointDVI()], [dviData()]
#' @export
#'
#' @examples
#' ### Family with three missing
#' familiasTodvir(famfile = "https://familias.name/dviapp/example8.fam")
#'
#' \dontrun{
#' ### No data for missing, an error is returned:
#' familiasTodvir("https://familias.name/dviapp/BrotherPower.fam")
#' }
#'
#' @importFrom forrel readFam
familiasTodvir = function(famfile, missingPrefix = "M", verbose = FALSE){
  
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
  
  if (length(x) == 2){ #cannot use lapply
    am = x[[2]][[2]]
    if(!is.ped(am)) am = am[[which(pedsize(am) > 1)]] # remove redundant singletons
  }
  else if (length(x) > 2) {
    am = lapply(x[-1], function(dat) {
      ref = dat$`Reference pedigree`
      ref = dat[[2]]
      if(!is.ped(ref))
        ref = ref[[which(pedsize(ref) > 1)]]     
      ref
    })
  }
  
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
  if(is.ped(am))
    dvi = relabelDVI(dvi0, missingPrefix = missingPrefix)
  else
    dvi = relabelDVI(dvi0, missingFormat = "M.[FAM].[IDX]")
  
  dvi
  
  # # Treat cases with one reference family and several separately
  # if (is.ped(am)){
  #   missing = grep("^Missing", unlist(labels(am)), value = TRUE)
  #   m = length(missing)
  #   if(m == 0)
  #     stop2("Reference pedigree without missing.")
  #   else{
  #     newMissing = paste(missingPrefix, 1:m, sep = ".")
  #     am = relabel(am, new = newMissing, old =  missing)
  #   }
  #   missing = newMissing
  # }
  # else if (is.pedList(am)){
  #   missing = NULL
  #   j = 1
  #   for (i in 1:length(am)){
  #     missing1 = grep("^Missing", unlist(labels(am[[i]])), value = TRUE)
  #     m = length(missing1)
  #     if(m == 0)
  #       stop2("Reference pedigree without missing.")
  #     else{
  #       for (k in 1:m){
  #         newMissing = paste(missingPrefix,j, sep = ".",k)
  #         am[[i]] = relabel(am[[i]], newMissing, missing1[k])
  #         missing = c(missing, newMissing)
  #       }
  #       j = j + m
  #     }
  #   }
  # }
  #   dviData(pm = pm, am = am, missing = missing)
}