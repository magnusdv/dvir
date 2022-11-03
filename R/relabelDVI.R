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
#' @param pm A list of singletons.
#' @param am A list of pedigrees.
#' @param missing Character vector with names of missing persons.
#' @param victimPrefix Prefix used to label PM individuals.
#' @param familyPrefix Prefix used to label the AM families.
#' @param refPrefix Prefix used to label the reference individuals, i.e., the
#'   typed members of the AM families.
#' @param missingPrefix Prefix used to label the missing persons in the AM
#'   families. The word "family" is treated as a special case, where the family
#'   name is used as prefix in each family, e.g., F1-1, F1-2, F2-1, ...
#' @param othersPrefix Prefix used to label other untyped individuals. Default:
#'   1, 2, ...
#'
#' @return A list with entries "pm", "am" and "missing".
#'
#' @examples
#'
#' # Builtin dataset `example2`
#' pm = example2$pm
#' am = example2$am
#' missing = example2$missing
#'
#' relabelDVI(pm, am, missing,
#'            victimPrefix  = "vic",
#'            familyPrefix  = "fam",
#'            refPrefix     = "ref",
#'            missingPrefix = "mp")
#'
#' # Family-wise numbering of missing persons
#' relabelDVI(pm, am, missing, missingPrefix = "family")
#'
#' @export
relabelDVI = function(pm, am, missing, 
                        victimPrefix = "V", familyPrefix = "F",
                        refPrefix = "R", missingPrefix = "M", othersPrefix = "") {
  
  if(is.singleton(pm)) 
    pm = list(pm)
  if(is.ped(am)) 
    am = list(am)
  
  npm = length(pm)
  nam = length(am)
  nmiss = length(missing)
  
  
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
  if(!is.null(missingPrefix)) {
  
    if(length(missingPrefix) != 1)
      stop2("`missingPrefix` must have length 1: ", missingPrefix)
    
    # Preliminary fix if same label is used in multiple components (e.g. "Missing person")
    if(nmiss == 1) {
      ocs = sapply(labels(am), function(v) missing %in% v)
      if(sum(ocs) > 1) {
        for(i in which(ocs)) 
          am[[i]] = relabel(am[[i]], missing, new = paste(missing, i))
        missing = paste(missing, which(ocs))
        nmiss = length(missing)
      }
    }
        
    # Special case
    if(identical(missingPrefix, "family")) {
      newmiss = character(nmiss)
      names(newmiss) = missing
      
      for(famname in names(am)) {
        m = intersect(missing, labels(am[[famname]]))
        newmiss[m] = paste(famname, seq_along(m), sep = "-")
      }
    } 
    else {
      # Default case
      newmiss = paste0(missingPrefix, 1:nmiss)
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
  
  
  # Return new objects
  list(pm = pm, am = am, missing = missing)
}

