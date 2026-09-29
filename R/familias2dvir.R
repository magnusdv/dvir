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
                         verbose = TRUE, missingIdentifier = "^Missing"){
  
  # Read fam file
  x = pedFamilias::readFam(famfile, useDVI = TRUE, verbose = verbose)

  # PM data -----------------------------------------------------------------
  
  pm = x$`Unidentified persons`
  if(is.null(pm))
    stop2("No `Unidentified persons` found")

  # AM data -----------------------------------------------------------------
  
  am = x[-1]
  
  hasAM = length(am)
  inspectAM = character()
  
  if(hasAM) {
   
    # Remove untyped components; catch disconnected fams
    z = lapply(am, function(ref) {
      if(is.ped(ref))
         return(list(ped = ref, inspect = FALSE))
      
      # Remove Reference pedigree if present
      idx = match("Reference pedigree", names(ref), nomatch = 0)
      if(length(ref) == 2 && idx > 0)
        ref = ref[[3 - idx]] # the other
      
      # Expect single component with typed references
      cmp = getComponent(ref, typedMembers(ref))
      if(max(cmp) > min(cmp))
        list(ped = .connectPed(ref[unique.default(cmp)]), inspect = TRUE)
      else
        list(ped = ref[[cmp[1]]], inspect = FALSE)
    })
  
    inspectAM = names(z)[vapply(z, `[[`, logical(1), "inspect")]
    am = lapply(z, `[[`, "ped")
    
    # Check for duplicated names among reference individuals
    if(length(am) > 0) {
      typed = typedMembers(am)
      dups = anyDuplicated.default(typed)
      if(dups) {
        if(verbose) cat(sprintf("Warning: Duplicated reference name '%s'\n", typed[dups]), 
                        "        Adding family prefix\n")
        for(nam in names(am)) {
          y = am[[nam]]
          old = intersect(typed, y$ID)
          am[[nam]] = relabel(y, old = old, new = paste(nam, old, sep = "_"))
        }
      }
    }
  }
  
  # Missing individuals -----------------------------------------------------

  if(hasAM) {
    missing = grep(missingIdentifier, unlist(labels(am)), value = TRUE)
    
    if(!length(missing) && verbose)
      cat("No missing individuals indicated\n")
  }  
  else missing = NULL

  # Create DVI object -------------------------------------------------------

  dvi0 = dviData(pm = pm, am = am, missing = missing, generatePairings = FALSE)
  
  # Relabel missing persons (NB: pairings are generated here)
  dvi = relabelDVI(dvi0, victimPrefix = victimPrefix, familyPrefix = familyPrefix,
             refPrefix = refPrefix, missingPrefix = missingPrefix, 
             missingFormat = missingFormat, othersPrefix = othersPrefix)
  
  if(length(inspectAM))
    attr(dvi, "reconnectedAM") = inspectAM
  
  dvi
}



# Utility for "fixing" disconnected reference families
# Typically these appear as [main comp] + [spouse]
.connectPed = function(x) {
  if(!is.pedList(x) || length(x) != 2L)
    stop2("Input must contain exactly two pedigree components")

  if(is.singleton(x[[1]]))
    x = x[2:1]
  
  x1 = x[[1]]
  x2 = x[[2]]
  if("Missing person" %notin% x1$ID)
    stop2("Expected a pedigree member named 'Missing person'")

  p = c("Missing person", x2$ID)
  sx = getSex(x, p)
  if(sx[1] != sx[2])
     return(addChild(x, p, id = "dummy", sex = 0, verbose = FALSE))
  
  # Same sex: connect through a shared spouse
  x |> 
    addChild(c(p[1], "_dummy1"), id = "_dummy2", sex = 0, verbose = FALSE) |> 
    addChild(c("_dummy1", p[2]), id = "_dummy3", sex = 0, verbose = FALSE)
}