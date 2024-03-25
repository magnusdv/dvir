#' Automatic labelling of a DVI dataset
#'
#' Relabel the individuals and families in a DVI dataset.
#'
#' @param dvi A `dviData` object, typically created with [dviData()].
#' @param victims A named vector of the form `c(old = new)` with names for the
#'   PM samples, or a function to be applied to the existing names.
#' @param victimPrefix Prefix used to label PM individuals.
#' @param familyPrefix Prefix used to label the AM families.
#' @param refs A named vector of the form `c(old = new)` with names for the
#'   typed references, or a function to be applied to the existing names.
#' @param refPrefix Prefix used to label the reference individuals, i.e., the
#'   typed members of the AM families.
#' @param missingPrefix Prefix used to label the missing persons. At most one of
#'   `missingPrefix` and `missingFormat` can be given.
#' @param missingFormat A string indicating family-wise labelling of missing
#'   persons, using `[FAM]`, `[IDX]`, `[MIS]` as place holders with the
#'   following meanings (see Examples):
#'   * `[FAM]`: family index
#'   * `[IDX]`: index of missing person within the family
#'   * `[MIS]`: index within all missing persons
#' @param othersPrefix Prefix used to label other untyped individuals. Use ""
#'   for numeric labels ( 1, 2, ...).
#'
#' @return A [dviData()] object.
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
relabelDVI = function(dvi, victims = NULL, victimPrefix = NULL, familyPrefix = NULL,
                      refs = NULL, refPrefix = NULL, missingPrefix = NULL, 
                      missingFormat = NULL, othersPrefix = NULL) {
  
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
    
    victims = paste0(victimPrefix, 1:npm)
  }
  
  if(!is.null(victims)) {
    pm = relabel(pm, new = victims)
  }
  
  # Relabel typed reference individuals
  if(!is.null(refPrefix)) {
    
    if(length(refPrefix) != 1)
      stop2("`refPrefix` must have length 1: ", refPrefix)
    
    refs = paste0(refPrefix, seq_along(typedMembers(am)))
  }
  
  if((!is.null(refs))) {
    am = relabel(am, new = refs, old = typedMembers(am))
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
    
    # Temporary fix if same label is used in multiple components (e.g. "Missing person")
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
        stop2("`missingFormat` must have length 1: ", missingFormat)
      
      mis = seq_len(nmiss)
      fam = getComponent(am, missing)
      idx = integer(length(missing))
      for(i in unique.default(fam))
        idx[fam == i] = seq_len(sum(fam == i))
      
      fmt = sub("[FAM]", "%{fam}d", fixed = TRUE,
                sub("[IDX]", "%{idx}d", fixed = TRUE,
                    sub("[MIS]", "%{mis}d", missingFormat, fixed = TRUE)))
      newmiss = sprintfNamed(fmt, fam = fam, idx = idx, mis = mis)
      
      if(length(newmiss) != nmiss)
        newmiss = rep_len(newmiss, length.out = nmiss)
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

