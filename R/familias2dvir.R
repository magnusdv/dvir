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
   
    # Repair disconnected reference families
    z = lapply(am, function(ref) {
      if(is.ped(ref))
        return(list(ped = ref, inspect = character()))

      # Remove Reference pedigree if present
      idx = match("Reference pedigree", names(ref), nomatch = 0)
      if(length(ref) == 2 && idx > 0)
        ref = ref[[3 - idx]] # the other
      if(is.ped(ref))
        return(list(ped = ref, inspect = character()))

      # Keep only components with references or missing persons
      ref = Filter(function(p) length(typedMembers(p)) ||
                     any(grepl(missingIdentifier, p$ID)), ref)
      if(!length(ref))
        return(NULL)
      if(length(ref) == 1L)
        return(list(ped = ref[[1]], inspect = character()))

      # Put component containing Missing person first
      mp = getComponent(ref, "Missing person", checkUnique = TRUE)
      if(!is.na(mp))
        ref = ref[c(mp, seq_along(ref)[-mp])]

      mpSex = getSex(ref[[1]], "Missing person")
      ref2 = .repairSpouse(ref)
      nonSpouse = !is.ped(ref2)
      
      # If still not connected, try to connect Child components
      if(nonSpouse)
        ref2 = .repairChild(ref2)
        
      # If still not connected, resort to ad hoc connection
      if(!is.ped(ref2))
        ref2 = .connectPed(ref2)
      
      list(ped = ref2, inspect = c(
        if(mpSex %in% 1:2 && getSex(ref2, "Missing person") != mpSex) "mpSexChanged",
        if(nonSpouse) "nonSpouseRepair"))
    })
    z = Filter(Negate(is.null), z)

    inspectAM = lapply(z, `[[`, "inspect")
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
    if(anyDuplicated(missing) && is.null(missingPrefix) && is.null(missingFormat))
      missingPrefix = "M"

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
  
  if(any(lengths(inspectAM))) {
    reasons = c("mpSexChanged", "nonSpouseRepair")
    attr(dvi, "inspectAM") = setNames(lapply(reasons, function(r)
      names(dvi$am)[vapply(inspectAM, function(x) r %in% x, logical(1))]), reasons)
  }
  
  dvi
}


# Utility for repairing disconnected spouse components
.repairSpouse = function(x) {
  if("Missing person" %notin% x[[1]]$ID)
    return(x)

  # Find explicit spouse in another component
  sp = labels(x)
  sp = sp[grepl("spouse|wife", sp, ignore.case = TRUE)]
  if(length(sp) != 1L || (spComp = getComponent(x, sp)) == 1L)
    return(x)

    # Find artificial spouses on both sides
  mpCopy = spouses(x[[spComp]], sp)
  mpCopy = mpCopy[startsWith(mpCopy, "added_")]
  spCopy = spouses(x[[1]], "Missing person")
  spCopy = spCopy[startsWith(spCopy, "added_")]
  if(!length(mpCopy) && !length(spCopy))
    return(x)

  # Correct sex of Missing person if needed
  if(length(mpCopy)) {
    mpSex = unique.default(getSex(x[[spComp]], mpCopy))
    if(length(mpSex) != 1L || mpSex == 0)
      return(x)
    if(getSex(x[[1]], "Missing person") == 0)
      x[[1]] = setSex(x[[1]], "Missing person", mpSex)
    else if(getSex(x[[1]], "Missing person") != mpSex)
      x[[1]] = swapSex(x[[1]], "Missing person", verbose = FALSE)
  }

  # Correct sex of Spouse if needed
  if(length(spCopy)) {
    spSex = unique.default(getSex(x[[1]], spCopy))
    if(length(spSex) != 1L || spSex == 0)
      return(x)
    if(getSex(x[[spComp]], sp) == 0)
      x[[spComp]] = setSex(x[[spComp]], sp, spSex)
    else if(getSex(x[[spComp]], sp) != spSex)
      x[[spComp]] = swapSex(x[[spComp]], sp, verbose = FALSE)
  }

  # Replace artificial spouses by their real counterparts
  .replaceParents(x, c(mpCopy, spCopy),
                  c(rep("Missing person", length(mpCopy)), rep(sp, length(spCopy))))
}

# Replace artificial parents of Child* with Missing person
.repairChild = function(x) {
  if(length(x) != 2L || "Missing person" %notin% x[[1]]$ID)
    return(x)

  y = x[[2]]
  ch = grep("^Child", y$ID, value = TRUE)
  if(!length(ch))
    return(x)
  copies = lapply(ch, function(id) grep("^added_", parents(y, id), value = TRUE))
  if(any(lengths(copies) != 1L))
    return(x)
  copies = unique.default(unlist(copies, use.names = FALSE))
  sx = unique.default(getSex(y, copies))
  if(length(sx) != 1L || sx == 0)
    return(x)
  mpSex = getSex(x[[1]], "Missing person")
  if(mpSex == 0)
    x[[1]] = setSex(x[[1]], "Missing person", sx)
  else if(mpSex != sx)
    x[[1]] = swapSex(x[[1]], "Missing person", verbose = FALSE)

  .replaceParents(x, copies, rep("Missing person", length(copies)))
}

# Replace artificial parents and remove them
.replaceParents = function(x, old, new) {
  m = as.matrix(x)
  labs = attr(m, "LABELS")
  par = m[,2:3, drop = FALSE]
  idx = match(par, match(old, labs))
  par[!is.na(idx)] = match(new, labs)[idx[!is.na(idx)]]
  m[,2:3] = par
  keep = labs %notin% old
  attrs = attributes(m)
  attrs$LABELS = labs[keep]
  restorePed(m[keep, , drop = FALSE], attrs = attrs)
}

# Link disconnected components with untyped children
.connectPed = function(x) {
  Reduce(function(a, b) {
    p = c(if("Missing person" %in% a$ID) "Missing person" else a$ID[1], b$ID[1])
    if(getSex(a, p[1]) == 0)
      a = setSex(a, p[1], 2)
    if(getSex(b, p[2]) == 0)
      b = setSex(b, p[2], 2)
    if(getSex(a, p[1]) == getSex(b, p[2])) {
      a = addChild(a, c(p[1], "_dummy"), sex = 0, verbose = FALSE)
      p[1] = "_dummy"
    }
    addChild(list(a, b), p, sex = 0, verbose = FALSE)
  }, x)
}
