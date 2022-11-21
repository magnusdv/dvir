#' Automatic labelling of a DVI dataset
#'
#' Relabel the families and individuals in a DVI dataset, using automatic
#' labelling.
#'
#' By default, the following labelling scheme is applied:
#'
#' * Victims (PM data): V1, V2, ...
#'
#' * Reference families: F1, F2, ...
#'
#' * Reference individuals: R1, R2, ...
#'
#' * Missing persons: M1, M2, ...
#'
#' * Others: 1, 2, ...
#'
#' @param dvi A `dviData` object, typically created with [dviData()].
#' @param victimPrefix Prefix used to label PM individuals.
#' @param familyPrefix Prefix used to label the AM families.
#' @param refPrefix Prefix used to label the reference individuals, i.e., the
#'   typed members of the AM families.
#' @param missingPrefix Prefix used to label the missing persons in the AM
#'   families. The word "family" is treated as a special case, where the family
#'   name is used as prefix in each family, e.g., F1-1, F1-2, F2-1, ...
#' @param missingFormat A string indicating family-wise labelling of missing
#'   persons, using `[FAM]` an `[IDX]` as place holders for the family index and
#'   the missing person index within the family. See Examples.
#' @param othersPrefix Prefix used to label other untyped individuals. Default:
#'   1, 2, ...
#'
#' @return A list with entries "pm", "am" and "missing".
#'
#' @examples
#'
#' # Builtin dataset `example2`
#' relabelDVI(example2,
#'            victimPrefix  = "vic",
#'            familyPrefix  = "fam",
#'            refPrefix     = "ref",
#'            missingPrefix = "mp")
#'
#' # Family-wise labelling of missing persons
#' relabelDVI(example2, missingFormat = "M[FAM]-[IDX]")
#' relabelDVI(example2, missingFormat = "M[IDX] (F[FAM])")
#' relabelDVI(example2, missingFormat = "fam[FAM].m[IDX]")
#'
#' @export
relabelDVI = function(dvi, victimPrefix = NULL, familyPrefix = NULL,
                      refPrefix = NULL, missingPrefix = NULL, 
                      missingFormat = NULL,
                      othersPrefix = NULL) {
  
  dvi = consolidateDVI(dvi)
  
  pm = dvi$pm
  am = dvi$am
  missing = dvi$missing
  
  npm = length(pm)
  nam = length(am)
  nmiss = length(missing)
  
  if(identical(missingPrefix, "family") && is.null(names(am)) && is.null(familyPrefix))
    familyPrefix = "F"
  
  # Relabel PM data
  if(!is.null(victimPrefix)) {
    
    if(length(victimPrefix) != 1)
      stop2("`victimPrefix` must have length 1: ", victimPrefix)
    
    newvics = paste0(victimPrefix, 1:npm)
    pm = relabel(pm, newvics)
    names(pm) = newvics
  }
  
  # Relabel typed reference individuals
  if(!is.null(refPrefix)) {
    
    if(length(refPrefix) != 1)
      stop2("`refPrefix` must have length 1: ", refPrefix)
    
    refs = typedMembers(am)
    nrefs = length(refs)
    am = relabel(am, old = refs, new = paste0(refPrefix, 1:nrefs))
  }
  
  # Rename AM families
  if(!is.null(familyPrefix)) {
    
    if(length(familyPrefix) != 1)
      stop2("`familyPrefix` must have length 1: ", familyPrefix)
    
    names(am) = paste0(familyPrefix, 1:nam)
  }
  
  
  # Relabel missing persons
  if(!is.null(missingPrefix) || !is.null(missingFormat)) {
  
    if(!is.null(missingPrefix) && !is.null(missingFormat))
      stop2("At most one of `missingPrefix` and `missingFormat` can be used")
    
    # Preliminary fix if same label is used in multiple components (e.g. "Missing person")
    labs = unlist(labels(am), use.names = FALSE)
    if(any(labs[duplicated(labs)] %in% missing)) {
      tmpmiss = character()
      
      for(i in seq_along(am)) {
        comp = am[[i]]
        missi = intersect(missing, comp$ID)
        if(!length(missi))
          next
        am[[i]] = relabel(comp, missi, new = paste(missi, i, sep = "."))
        tmpmiss = c(tmpmiss, paste(missi, i, sep = "."))
      }
      
      missing = tmpmiss
      nmiss = length(missing)
    }
    
    # Regular prefix format
    if(!is.null(missingPrefix)){
      
      if(length(missingPrefix) != 1)
        stop2("`missingPrefix` must have length 1: ", missingPrefix)
      
      # Default case
      newmiss = paste0(missingPrefix, 1:nmiss)
    }    
    
    # Special format
    if(!is.null(missingFormat)) {
      
      if(length(missingFormat) != 1)
        stop2("`missingPrefix` must have length 1: ", missingPrefix)
      
      fam = getComponent(am, missing)
      idx = integer(length(missing))
      for(i in unique.default(fam))
        idx[fam == i] = seq_len(sum(fam == i))
      
      # Location of [FAM] and [IDX] in the string
      famLoc = regexpr("[FAM]", missingFormat, fixed = T)
      idxLoc = regexpr("[IDX]", missingFormat, fixed = T)
      if(famLoc == -1 || idxLoc == -1)
        stop2("`missingFormat` should be a string containing '[FAM]' and '[IDX]'")
      
      format = sub("[FAM]", "%d", sub("[IDX]", "%d", missingFormat, fixed = TRUE), fixed = TRUE)
      if(famLoc < idxLoc)
        newmiss = sprintf(format, fam, idx)
      else
        newmiss = sprintf(format, idx, fam)
    } 
    
    # Relabel
    am = relabel(am, old = missing, new = newmiss)
    
    # Update missing vector
    missing = as.character(newmiss)
  }
  
  # Other untyped individuals
  if(!is.null(othersPrefix)) {
    
    if(length(othersPrefix) != 1)
      stop2("`othersPrefix` must have length 1: ", othersPrefix)
  
    k = 0
    for(i in 1:nam) {
      oth = setdiff(untypedMembers(am[[i]]), missing)
      if(length(oth)) {
        am[[i]] = relabel(am[[i]], old = oth, new = paste0(othersPrefix, k + seq_along(oth)))
        k = k + length(oth)
      }
    }
  }
  
  # Return new object
  dviData(pm = pm, am = am, missing = missing)
}

