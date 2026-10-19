#' Automatic labelling of a DVI dataset
#'
#' Relabel the individuals and families in a DVI dataset.
#'
#' @param dvi A `dviData` object, typically created with [dviData()].
#' @param victims A named vector `c(old = "new")` for PM samples, or a function applied to
#'   the existing labels.
#' @param victimPrefix Prefix used to label PM individuals.
#' @param families New family names (in order), a named vector `c(old = "new")`, or a
#'   function applied to the existing family names. Overrides `familyPrefix`.
#' @param familyPrefix Prefix used to label AM families.
#' @param refs A named vector `c(old = "new")` for typed references, or a function applied
#'   to their existing labels.
#' @param refPrefix Prefix used to label typed reference individuals. Ignored if `refs` is
#'   given.
#' @param missingPrefix Prefix used to label missing persons. Cannot be combined with
#'   `missingFormat`.
#' @param missingFormat Format for missing-person labels; see Details.
#' @param othersPrefix Prefix used to label other untyped individuals. Use "" for numeric
#'   labels (1, 2, ...).
#'
#' @details `missingFormat` supports these keywords:
#'
#' * `[FIDX]`: family index in the AM list
#' * `[FNAME]`: family name
#' * `[MIDX]`: index among all missing persons
#' * `[MIDX_FAM]`: index of missing person within its family
#'
#' The following deprecated keywords are still supported, for back compatibility:
#' 
#' * `[FAM]` (same as `[FIDX]`)
#' * `[MIS]` (same as `[MIDX]`)
#' * `[IDX]` (same as `[MIDX_FAM]`)
#'
#' @return A [dviData()] object.
#'
#' @examples
#' relabelDVI(example2, victimPrefix = "vic", familyPrefix = "fam",
#'            refPrefix = "ref", missingPrefix = "mp")
#'
#' relabelDVI(example2, missingFormat = "Missing[MIDX]")
#' relabelDVI(example2, missingFormat = "M[FIDX]-[MIDX_FAM]")
#' relabelDVI(example2, familyPrefix = "Case", missingFormat = "[FNAME]-M[MIDX_FAM]")
#'
#' @export
relabelDVI = function(dvi, victims = NULL, victimPrefix = NULL, familyPrefix = NULL,
                      refs = NULL, refPrefix = NULL, missingPrefix = NULL,
                      missingFormat = NULL, othersPrefix = NULL, families = NULL) {
  
  dvi = consolidateDVI(dvi)
  
  pm = dvi$pm
  am = dvi$am
  missing = dvi$missing
  
  npm = length(pm)
  nam = length(am)
  
  # Relabel PM data
  if(npm > 0) {
    if(!is.null(victims))
      pm = relabel(pm, new = victims)
    else if(!is.null(victimPrefix)) {
      if(length(victimPrefix) != 1)
        stop2("`victimPrefix` must have length 1: ", victimPrefix)
      victims = paste0(victimPrefix, 1:npm)
      pm = relabel(pm, new = victims)
    }
  }
  
  # If no AM data: Return early
  if(nam == 0)
    return(dviData(pm = pm, am = am, missing = missing))  
  

  # AM data -------------------------------------------------------------------------------------

  # Relabel typed reference individuals
  if(is.null(refs) && !is.null(refPrefix)) {
    if(length(refPrefix) != 1)
      stop2("`refPrefix` must have length 1: ", refPrefix)
    refs = paste0(refPrefix, seq_along(typedMembers(am)))
  }
  if(!is.null(refs))
    am = relabel(am, new = refs,
                 old = if(is.null(names(refs))) typedMembers(am) else names(refs))

  # Rename AM families
  if(!is.null(families)) {
    famnames = names(am)
    if(is.function(families))
      families = families(famnames)
    if(!is.null(names(families)))
      famnames[match(names(families), famnames)] = unname(families)
    else
      famnames = families
    if(length(famnames) != nam || anyNA(famnames) || anyDuplicated.default(famnames))
      stop2("Invalid or duplicated family names")
    names(am) = famnames
  }
  else if(!is.null(familyPrefix)) {
    if(length(familyPrefix) != 1)
      stop2("`familyPrefix` must have length 1: ", familyPrefix)
    names(am) = paste0(familyPrefix, seq_len(nam))
  }

  # Relabel missing persons
  if(length(missing) && (!is.null(missingPrefix) || !is.null(missingFormat))) {
    if(!is.null(missingPrefix) && !is.null(missingFormat))
      stop2("At most one of `missingPrefix` and `missingFormat` can be used")

    # Disambiguate missing-person labels repeated across families
    labs = unlist(labels(am), use.names = FALSE)
    if(any(labs[duplicated.default(labs)] %in% missing)) {
      tmpmiss = vector("list", nam)
      for(i in seq_along(am)) {
        missi = intersect(missing, am[[i]]$ID)
        if(length(missi)) {
          tmpmiss[[i]] = paste0(missi, ".", i)
          am[[i]] = relabel(am[[i]], old = missi, new = tmpmiss[[i]])
        }
      }
      missing = unlist(tmpmiss, use.names = FALSE)
    }

    # Regular prefix format
    if(!is.null(missingPrefix)) {
      if(length(missingPrefix) != 1)
        stop2("`missingPrefix` must have length 1: ", missingPrefix)
      newmiss = paste0(missingPrefix, seq_along(missing))
    }
    else {
      if(length(missingFormat) != 1)
        stop2("`missingFormat` must have length 1: ", missingFormat)

      midx = seq_along(missing)
      fidx = getComponent(am, missing)
      midx_fam = integer(length(missing))
      for(i in unique.default(fidx))
        midx_fam[fidx == i] = seq_len(sum(fidx == i))

      keys = c(FIDX = "%{fidx}d", FNAME = "%{fname}s",
               MIDX = "%{midx}d", MIDX_FAM = "%{midx_fam}d",
               FAM = "%{fidx}d", MIS = "%{midx}d", IDX = "%{midx_fam}d")
      fmt = missingFormat
      for(k in names(keys))
        fmt = gsub(paste0("[", k, "]"), keys[[k]], fmt, fixed = TRUE)
      
      bad = regmatches(fmt, regexpr("\\[[^]]+\\]", fmt))
      if(length(bad))
        stop2("Unknown keyword: ", substring(bad, 2L, nchar(bad) - 1L))
      newmiss = sprintfNamed(fmt, fidx = fidx, fname = names(am)[fidx],
                             midx = midx, midx_fam = midx_fam)
    }

    am = relabel(am, old = missing, new = newmiss)
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

