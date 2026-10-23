#' Construct AM pedigrees from a relationship table
#'
#' Each relationship is relative to one missing person per family. Marker data are
#' returned separately by default. Unspecified relatives are placed on the simplest shared
#' branch.
#'
#' @param x A data frame or path to a tab-separated file. The first three columns must be
#'   Family, Relationship (or Kinship), and Sample, in that order.
#' @param mpFormat Naming of missing persons: "M" gives M1, M2, ...; "family" gives F1_MP,
#'   F2_MP, ...
#' @param attachGeno If TRUE, attach genotypes to the pedigrees.
#' @param loci Optional list of marker attributes passed on to [setMarkersFromMatrix()].
#'   If NULL, observed alleles are assigned equal frequencies when attaching markers.
#'
#' @return A list with `am` (one pedigree or list of components per family), `missing` (MP
#'   labels), `issues` (family IDs grouped by issue), and `geno` (genotypes, one column
#'   per marker). If `attachGeno = TRUE`, the genotype data is attached to the `am`
#'   pedigrees. contains markers. Disconnected references remain singletons; construction
#'   errors return all singletons.
#'
#' @examples
#' x = data.frame(
#'   Family = rep(c("F1", "F2", "F3"), 1:3),
#'   Relationship = c("[Mother]", "[Wife]", "[Daughter]",
#'                    "[Grandmother]", "[Brother]", "[Niece]"),
#'   Sample = paste0("S", 1:6),
#'   M1 = c("1/1", "1/2", "2/2", "1/1", "1/2", "2/2"))
#' x
#' 
#' y = buildAM(x, attachGeno = TRUE)
#' 
#' plotPedList(y$am, hatched = typedMembers, marker = 1)
#' y$issues
#'
#' @export
#' @importFrom utils read.delim
buildAM = function(x, mpFormat = c("M", "family"), attachGeno = FALSE, loci = NULL) {
  mpFormat = match.arg(mpFormat)
  if(is.character(x) && length(x) == 1L)
    x = read.delim(x, check.names = FALSE, colClasses = "character")

  # Remove trailing empty columns without headers
  blank = !nzchar(names(x)) & vapply(x, \(z) all(is.na(z) | !nzchar(trimws(z))), TRUE)
  x = x[, seq_len(max(c(3L, which(!blank)))), drop = FALSE]
  
  a = list2DF(list(family = trimws(as.character(x[[1]])),
                   relation = trimws(as.character(x[[2]])),
                   sample = trimws(as.character(x[[3]]))))
  
  if(anyNA(a$relation) || any(!grepl("^\\[[^][]+\\]$", a$relation)))
    stop2("Relationships must contain bracketed keys")
  
  # Normalise relationships and make sample IDs unique across families
  a$relation = tolower(substring(a$relation, 2L, nchar(a$relation) - 1L))
  a$relation = trimws(gsub("[[:space:]]+", " ", gsub("[-_]", " ", a$relation)))
  
  # Add family prefixes unless all samples already have them
  haveFams = startsWith(a$sample, a$family)
  a$label = if(all(haveFams)) a$sample else paste0(a$family, "_", a$sample)

  dup = anyDuplicated.default(a$label)
  if(dup)
    stop2("Duplicated sample ID: ", a$label[dup])
  a$label = paste0(a$label, " [", a$relation, "]")

  # Collapse paired allele columns; keep combined genotype columns unchanged
  geno = x[, -(1:3), drop = FALSE]
  pairs = grep("[-_. ][12]$", names(geno))
  if(length(pairs)) {
    if(length(pairs) %% 2L)
      stop2("Odd number of allele columns")
    odd = pairs[c(TRUE, FALSE)]
    even = pairs[c(FALSE, TRUE)]
    markers = sub("[-_. ]1$", "", names(geno)[odd])
    if(!identical(markers, sub("[-_. ]2$", "", names(geno)[even])))
      stop2("Allele columns must have matching names, with 1/2 suffixes")
    vals = Map(\(u, v) paste(geno[[u]], geno[[v]], sep = "/"), odd, even)
    geno = geno[, -pairs, drop = FALSE]
    geno[markers] = vals
  }
  rownames(geno) = a$label

  families = unique.default(a$family)
  groups = split(a, factor(a$family, levels = families))
  mp = if(mpFormat == "M") paste0("M", seq_along(groups)) else paste0(families, "_MP")

  # Keep one pedigree or a list of components for each family
  z = lapply(seq_along(groups), function(i)
    tryCatch(.buildAMfamily(groups[[i]], mp[i]), error = function(e) 
      list(ped = singletons(c(mp[i], groups[[i]]$label), sex = 0),
           issues = conditionMessage(e))))
  
  names(z) = families
  am = lapply(z, `[[`, "ped")
  issue = lapply(z, `[[`, "issues")
  issues = split(rep(families, lengths(issue)), unlist(issue, use.names = FALSE))
  
  if(attachGeno && ncol(geno)) {
    
    # Split genotypes into two alleles per marker
    amat = splitGenotypes(geno, clean = TRUE, colNames = TRUE)
    
    # Match supplied locus attributes to genotype columns
    if(!is.null(loci)) {
      j = match(names(geno), vapply(loci, `[[`, "", "name"))
      if(anyNA(j)) 
        stop2("Missing locus attributes: ", names(geno)[is.na(j)])
      loci = loci[j]
    }
    am = setMarkersFromMatrix(am, amat, loci)
  }
  list(am = am, missing = unname(mp), issues = issues, geno = geno)
}


.buildAMfamily = function(a, mp) {
  family = a$family[1]
  r = a$relation
  # Classify relations by type, parental side and half-sibling status
  half = grepl("\\bhalf\\b", r, perl = TRUE)
  maternal = grepl("\\bmaternal\\b", r, perl = TRUE)
  paternal = grepl("\\bpaternal\\b", r, perl = TRUE)
  side = ifelse(maternal, "M", ifelse(paternal, "P", ""))
  type = trimws(gsub("[[:space:]]+", " ",
                      gsub("\\b(half|maternal|paternal)\\b", "", r, perl = TRUE)))
  sx = c(father = 1, mother = 2, parent = 0, son = 1, daughter = 2, child = 0,
         brother = 1, sister = 2, sibling = 0, grandfather = 1, grandmother = 2,
         grandson = 1, granddaughter = 2, grandchild = 0, uncle = 1, aunt = 2,
         niece = 2, nephew = 1, wife = 2, husband = 1, spouse = 0,
         direct = 0, identity = 0)
  valid = type %in% names(sx) & (!half | type %in% c("brother", "sister", "sibling")) &
    (side == "" | half | type %in% c("grandfather", "grandmother", "uncle", "aunt")) &
    !(maternal & paternal)

  # Keep unsupported and ignored references as separate singletons
  bad = !valid | type %in% c("direct", "identity")
  issues = c(if(any(!valid)) "Unsupported relationship",
             if(any(type %in% c("direct", "identity"))) "Ignored relationship")
  if(all(bad))
    return(list(ped = singletons(c(mp, a$label), sex = 0), issues = issues))
  
  ignored = a$label[bad]
  a = a[!bad, , drop = FALSE]
  r = a$relation
  type = type[!bad]
  half = half[!bad]
  side = side[!bad]
  refSex = unname(sx[type])
  n = length(type)

  sibling = type %in% c("brother", "sister", "sibling") & !half
  child = type %in% c("son", "daughter", "child")
  grandchild = type %in% c("grandson", "granddaughter", "grandchild")
  niece = type %in% c("niece", "nephew")
  uncle = type %in% c("uncle", "aunt")
  grandparent = type %in% c("grandfather", "grandmother")
  parent = type %in% c("father", "mother", "parent")
  ambiguous =
    (any(uncle & side == "") && (sum(uncle) > 1L || any(parent | grandparent | half))) ||
    (any(grandparent & side == "") && (sum(grandparent) > 1L || any(parent | uncle | half))) ||
    (any(half & side == "") && any(parent | uncle | grandparent)) ||
    (sum(half) > 1L && (any(side[half] == "") || anyDuplicated.default(side[half]) > 0L)) ||
    (any(type == "parent") && !any(type %in% c("father", "mother")) && any(uncle | grandparent | half)) ||
    (any(niece) && (sum(niece) > 1L || any(sibling | half))) ||
    (any(grandchild) && (sum(grandchild) > 1L || any(child)))
  
  if(ambiguous) 
    issues = c(issues, "Ambiguous pedigree")
  
  desc = any(child | grandchild)
  if(any(type %in% c("wife", "husband", "spouse")) && !desc)
    stop2("Spouse without children cannot be represented")
  
  mpSex = if(any(type == "wife")) 1L else if(any(type == "husband")) 2L else if(desc) 1L else 0L

  # Match references to structural positions before creating the pedigree
  slot = rep("", n)
  slot[type == "father"] = "F"
  slot[type == "mother"] = "M"
  slot[type %in% c("wife", "husband", "spouse")] = "SP"
  j = which(type == "parent")
  free = .mysetdiff(c("F", "M"), slot)
  if(length(j) > length(free)) stop2("Conflicting parents of MP")
  if(length(j)) slot[j[order(a$sample[j])]] = free[seq_along(j)]

  # Prefer the side of a known parent, otherwise the paternal side
  preferred = if(any(slot == "M") && !any(slot == "F")) "M" else "P"

  # Grandparents without a side favour an already represented branch
  j = which(grandparent)
  for(k in j[order(side[j] == "")]) {
    choices = if(type[k] == "grandfather") c("FF", "MF") else c("FM", "MM")
    if(side[k] == "P") choices = choices[1]
    if(side[k] == "M") choices = choices[2]
    if(side[k] == "" && (any(slot %in% c("MF", "MM")) && !any(slot %in% c("FF", "FM")) ||
                         !any(slot %in% c("FF", "FM", "MF", "MM")) && preferred == "M"))
      choices = rev(choices)
    free = .mysetdiff(choices, slot)
    if(!length(free)) 
      stop2("Incompatible grandparents")
    
    slot[k] = free[1]
  }

  # Half siblings share an outside parent; unqualified uncles share one side
  side[half & side == ""] = if(any(half & side == "M") && !any(half & side == "P")) "M" else preferred
  gpP = any(slot %in% c("FF", "FM"))
  gpM = any(slot %in% c("MF", "MM"))
  uSide = if(gpM && !gpP) "M" else if(gpP && !gpM) "P" else preferred
  side[uncle & side == ""] = uSide

  # A named sibling/child becomes parent of all nieces/grandchildren if available
  if(any(niece) && any(sibling | half)) {
    j = which(sibling | half)
    k = j[order(!sibling[j], match(type[j], c("brother", "sister", "sibling")), a$sample[j])][1]
    slot[k] = "SIB"
  }
  if(any(grandchild) && any(child)) {
    j = which(child)
    k = j[order(match(type[j], c("son", "daughter", "child")), a$sample[j])][1]
    slot[k] = "CH"
  }
  if(anyDuplicated.default(slot[nzchar(slot)])) 
    stop2("Conflicting reference relationships")

  # Only create ancestors and intermediate relatives actually needed
  anc = any(type %in% c("father", "mother", "parent", "brother", "sister", "sibling",
                        "grandfather", "grandmother", "uncle", "aunt", "niece", "nephew"))
  pat = any(slot %in% c("FF", "FM")) || any(uncle & side == "P")
  mat = any(slot %in% c("MF", "MM")) || any(uncle & side == "M")
  roles = c("MP", "F", "M", "FF", "FM", "MF", "MM", "SP", "SIB", "SIBSP", "CH", "CHSP", "PHM", "MHF")
  active = c(TRUE, anc, anc, pat, pat, mat, mat, desc, any(niece), any(niece),
             any(grandchild), any(grandchild), any(half & side == "P"), any(half & side == "M"))
  names(active) = roles
  id = setNames(paste0(family, "_", roles), roles)
  id["MP"] = mp
  id[slot[nzchar(slot)]] = a$label[nzchar(slot)]

  sibSex = if(any(slot == "SIB") && type[slot == "SIB"] == "sister") 2L else 1L
  chSex = if(any(slot == "CH") && type[slot == "CH"] == "daughter") 2L else 1L
  sex = c(MP = mpSex, F = 1L, M = 2L, FF = 1L, FM = 2L, MF = 1L, MM = 2L,
          SP = if(mpSex == 2) 1L else 2L, SIB = sibSex, SIBSP = 3L - sibSex,
          CH = chSex, CHSP = 3L - chSex, PHM = 2L, MHF = 1L)

  # Parent links for structural individuals
  fa = mo = setNames(rep("", length(roles)), roles)
  if(anc) {
    fa["MP"] = "F"
    mo["MP"] = "M"
  }
  if(pat) {
    fa["F"] = "FF"
    mo["F"] = "FM"
  }
  if(mat) {
    fa["M"] = "MF"
    mo["M"] = "MM"
  }
  if(any(niece)) {
    fromHalf = any(slot == "SIB") && half[slot == "SIB"]
    fa["SIB"] = if(fromHalf && side[slot == "SIB"] == "M") "MHF" else "F"
    mo["SIB"] = if(fromHalf && side[slot == "SIB"] == "P") "PHM" else "M"
  }
  if(any(grandchild)) {
    fa["CH"] = if(mpSex == 1) "MP" else "SP"
    mo["CH"] = if(mpSex == 1) "SP" else "MP"
  }

  # Parent links for references not occupying a structural position
  rf = rm = rep("", n)
  rf[sibling] = "F"
  rm[sibling] = "M"
  rf[half & side == "P"] = "F"
  rm[half & side == "P"] = "PHM"
  rf[half & side == "M"] = "MHF"
  rm[half & side == "M"] = "M"
  rf[uncle & side == "P"] = "FF"
  rm[uncle & side == "P"] = "FM"
  rf[uncle & side == "M"] = "MF"
  rm[uncle & side == "M"] = "MM"
  rf[child] = if(mpSex == 1) "MP" else "SP"
  rm[child] = if(mpSex == 1) "SP" else "MP"
  rf[niece] = if(sibSex == 1) "SIB" else "SIBSP"
  rm[niece] = if(sibSex == 1) "SIBSP" else "SIB"
  rf[grandchild] = if(chSex == 1) "CH" else "CHSP"
  rm[grandchild] = if(chSex == 1) "CHSP" else "CH"

  # Translate structural positions to person labels and construct once
  use = names(active)[active]
  extra = !nzchar(slot)
  p = ped(id = c(unname(id[use]), a$label[extra]),
          fid = unname(id[c(fa[use], rf[extra])]),
          mid = unname(id[c(mo[use], rm[extra])]),
          sex = c(unname(sex[use]), refSex[extra]), famid = family)
  if(!is.ped(p)) 
    stop2("Disconnected pedigree")
  
  added = .mysetdiff(p$ID, c(mp, a$label))
  if(length(added)) 
    p = relabel(p, old = added, new = as.character(seq_along(added)))
  if(length(ignored))
    p = c(list(p), singletons(ignored, sex = 0))
  
  list(ped = p, issues = issues)
}
