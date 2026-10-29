#' Construct AM pedigrees from a relationship table
#'
#' Each relationship is relative to one missing person per family. Marker data are
#' returned separately by default. Unspecified relatives are placed on the simplest shared
#' branch. Tags such as `[spouse:A]`, `[son:A]`, and `[paternal half brother:B]`
#' distinguish nuclear families and half-sibling groups. Untagged relatives
#' join a group if unique, otherwise remain singletons for manual review.
#'
#' @param x A data frame or path to a tab-separated file. The first three columns must be
#'   Family, Relationship (or Kinship), and Sample, in that order.
#' @param mpFormat Naming of missing persons:
#'   * "M" gives M1, M2, ...
#'   * "F" gives `<fam1>_MP`, `<fam2>_MP`, ..., where `<fam1>`, `<fam2>`, ... are the
#'   family names.
#' @param attachGeno If TRUE, attach genotypes to the pedigrees.
#' @param loci Optional list of marker attributes passed on to [setMarkersFromMatrix()].
#'   If NULL, observed alleles are assigned equal frequencies when attaching markers.
#'
#' @return A list with `am` (one pedigree or list of components per family), `missing` (MP
#'   labels), `issues` (family IDs grouped by issue), and `geno` (genotypes, one column
#'   per marker). If `attachGeno = TRUE`, the genotype data is attached to the `am`
#'   pedigrees. Disconnected references remain singletons; construction
#'   errors return all singletons.
#'
#' @examples
#' x = data.frame(
#'   Family = rep(c("F1", "F2", "F3"), 1:3),
#'   Relationship = c("[Mother]", "[Wife]", "[Daughter]",
#'                    "[Grandmother]", "[Brother]", "[Niece]"),
#'   Sample = paste0("S", 1:6),
#'   SNP1 = c("1/1", "1/2", "2/2", "1/1", "1/2", "2/2"))
#' x
#'
#' y = buildAM(x, attachGeno = TRUE)
#'
#' plotPedList(y$am, hatched = typedMembers, marker = 1)
#' y$issues
#'
#' # Two spouses with children from different nuclear families
#' z = data.frame(Family = "F4",
#'   Relationship = c("[spouse:A]", "[son:A]", "[spouse:B]", "[daughter:B]"),
#'   Sample = paste0("R", 1:4))
#' buildAM(z)$am
#'
#' # Half-niece through an unnamed half-sibling
#' w = data.frame(Family = "F5", Relationship = "[half niece]", Sample = "R5")
#' buildAM(w)$am
#'
#' @export
#' @importFrom utils read.delim
buildAM = function(x, mpFormat = c("M", "F"), attachGeno = FALSE, loci = NULL) {
  mpFormat = match.arg(mpFormat)
  if(is.character(x) && length(x) == 1L)
    x = read.delim(x, check.names = FALSE, colClasses = "character")

  # Remove trailing empty columns without headers
  blank = which(!nzchar(names(x)))
  if(length(blank)) {
    empty = vapply(x[blank], \(z) all(is.na(z) | !nzchar(trimws(z))), TRUE)
    x = x[, seq_len(max(c(3L, setdiff(seq_along(x), blank[empty])))), drop = FALSE]
  }
  
  a = list2DF(list(family = trimws(as.character(x[[1]])),
                   relation = trimws(as.character(x[[2]])),
                   sample = trimws(as.character(x[[3]]))))
  
  if(anyNA(a$relation) || any(!grepl("^\\[[^][]+\\]$", a$relation)))
    stop2("Relationships must contain bracketed keys")
  
  # Split nuclear-family tags before normalising relationship names
  a$relation = tolower(substring(a$relation, 2L, nchar(a$relation) - 1L))
  tagged = grepl(":", a$relation, fixed = TRUE)
  a$group = rep("", nrow(a))
  a$group[tagged] = toupper(trimws(sub("^[^:]*:", "", a$relation[tagged])))
  a$badGroup = tagged & !grepl("^[A-Z0-9]+$", a$group)
  a$relation = trimws(gsub("[[:space:]]+", " ",
                           gsub("[-_]", " ", sub(":.*$", "", a$relation))))

  # Parse all relationships once, rather than separately for each family
  a$half = grepl("\\bhalf\\b", a$relation, perl = TRUE)
  maternal = grepl("\\bmaternal\\b", a$relation, perl = TRUE)
  paternal = grepl("\\bpaternal\\b", a$relation, perl = TRUE)
  a$side = rep("", nrow(a))
  a$side[paternal] = "P"
  a$side[maternal] = "M"
  a$badSide = maternal & paternal
  a$type = trimws(gsub("[[:space:]]+", " ",
                       gsub("\\b(half|maternal|paternal)\\b", "", a$relation, perl = TRUE)))
  
  # Add family prefixes unless all samples already have them
  haveFams = startsWith(a$sample, a$family)
  a$label = if(all(haveFams)) a$sample else paste0(a$family, "_", a$sample)

  dup = anyDuplicated.default(a$label)
  if(dup)
    stop2("Duplicated sample ID: ", a$label[dup])
  suffix = rep("", nrow(a))
  suffix[tagged] = paste0(":", a$group[tagged])
  a$label = paste0(a$label, " [", a$relation, suffix, "]")

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
  mp = switch(mpFormat, M = paste0("M", seq_along(groups)), F = paste0(families, "_MP"))

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
  type = a$type
  half = a$half
  side = a$side
  sx = c(father = 1, mother = 2, parent = 0, son = 1, daughter = 2, child = 0,
         brother = 1, sister = 2, sibling = 0, grandfather = 1, grandmother = 2,
         grandson = 1, granddaughter = 2, grandchild = 0, uncle = 1, aunt = 2,
         niece = 2, nephew = 1, wife = 2, husband = 1, spouse = 0,
         direct = 0, identity = 0)
  spouse = type %in% c("wife", "husband", "spouse")
  child = type %in% c("son", "daughter", "child")
  halfSib = half & type %in% c("brother", "sister", "sibling")
  halfNiece = half & type %in% c("niece", "nephew")
  nuclear = spouse | child
  valid = type %in% names(sx) & (!half | halfSib | halfNiece) &
    (side == "" | half | type %in% c("grandfather", "grandmother", "uncle", "aunt")) &
    !a$badSide & !a$badGroup & (!nzchar(a$group) | nuclear | halfSib)

  # Match untagged relatives to a unique nuclear or half-sibling group
  bad = !valid | type %in% c("direct", "identity")
  ambiguousTags = rep(FALSE, length(type))
  for(s in c("", "P", "M")) {
    j = !bad & (if(s == "") nuclear else halfSib & side == s)
    grps = unique.default(a$group[j & nzchar(a$group)])
    if(length(grps) == 1L)
      a$group[j & !nzchar(a$group)] = grps
    if(length(grps) > 1L)
      ambiguousTags[j & !nzchar(a$group)] = TRUE
  }
  ambiguousGroup = ambiguousTags & nuclear
  ambiguousHalf = ambiguousTags & halfSib
  bad = bad | ambiguousTags

  # Half-nieces/nephews share a sibling, using an explicit side if available
  incompatibleNiece = rep(FALSE, length(type))
  j = which(!bad & halfNiece)
  if(length(j)) {
    known = side[j][nzchar(side[j])]
    nieceSide = if(length(known)) known[1] else ""
    incompatibleNiece = !bad & type %in% c("niece", "nephew") &
      (!halfNiece | (nzchar(side) & nzchar(nieceSide) & side != nieceSide))
    side[halfNiece & !bad & !nzchar(side)] = nieceSide
  }
  bad = bad | incompatibleNiece

  # Leave spouses without an identifiable child as singletons
  orphan = !bad & spouse & !(a$group %in% a$group[!bad & child])
  if(any(type %in% c("grandson", "granddaughter", "grandchild") & !bad) &&
     !any(child & !bad) && sum(spouse & !bad) == 1L)
    orphan[] = FALSE
  bad = bad | orphan
  issues = c(if(any(!valid)) "Unsupported relationship",
             if(any(type %in% c("direct", "identity"))) "Ignored relationship",
             if(any(ambiguousGroup)) "Ambiguous nuclear family",
             if(any(ambiguousHalf)) "Ambiguous half-sibling group",
             if(any(incompatibleNiece)) "Incompatible nieces",
             if(any(orphan)) "Spouse without children")
  if(all(bad))
    return(list(ped = singletons(c(mp, a$label), sex = 0), issues = issues))
  
  ignored = a$label[bad]
  a = a[!bad, , drop = FALSE]
  type = type[!bad]
  half = half[!bad]
  halfSib = halfSib[!bad]
  halfNiece = halfNiece[!bad]
  side = side[!bad]
  refSex = unname(sx[type])
  n = length(type)

  sibling = type %in% c("brother", "sister", "sibling") & !half
  spouse = type %in% c("wife", "husband", "spouse")
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
    (sum(halfSib) > 1L && (any(side[halfSib] == "") ||
                         anyDuplicated.default(side[halfSib & !nzchar(a$group)]) > 0L)) ||
    (any(type == "parent") && !any(type %in% c("father", "mother")) && any(uncle | grandparent | half)) ||
    (any(niece) && (sum(niece) > 1L || any(sibling | halfSib))) ||
    (any(grandchild) && (sum(grandchild) > 1L || any(child)))
  
  if(ambiguous) 
    issues = c(issues, "Ambiguous pedigree")
  
  desc = any(child | grandchild)
  if(any(type == "wife") && any(type == "husband"))
    stop2("Conflicting spouse sexes")
  mpSex = if(any(type == "wife")) 1L else if(any(type == "husband")) 2L else if(desc) 1L else 0L

  # Assign a structural spouse slot to each nuclear family
  spouseSlot = rep("SP", n)
  j = nzchar(a$group)
  spouseSlot[j] = paste0("SP_", a$group[j])
  spSlots = unique.default(spouseSlot[spouse | child])
  if(any(grandchild) && !length(spSlots)) spSlots = "SP"

  # Match references to structural positions before creating the pedigree
  slot = rep("", n)
  slot[type == "father"] = "F"
  slot[type == "mother"] = "M"
  slot[spouse] = spouseSlot[spouse]
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

  # Prefer a compatible named sibling/child as parent of nieces/grandchildren
  if(any(niece) && any(sibling | halfSib)) {
    j = which(sibling | halfSib)
    if(any(halfNiece))
      j = j[halfSib[j] & side[j] == side[which(halfNiece)[1]]]
    if(length(j)) {
      k = j[order(!sibling[j], match(type[j], c("brother", "sister", "sibling")), a$sample[j])][1]
      slot[k] = "SIB"
    }
  }
  if(any(grandchild) && any(child)) {
    j = which(child)
    k = j[order(match(type[j], c("son", "daughter", "child")), a$sample[j])][1]
    slot[k] = "CH"
  }
  if(anyDuplicated.default(slot[nzchar(slot)])) 
    stop2("Conflicting reference relationships")

  # Half-siblings in different groups have different outside parents
  halfSlot = ifelse(side == "P", "PHM", "MHF")
  j = halfSib & nzchar(a$group)
  halfSlot[j] = paste0(halfSlot[j], "_", a$group[j])
  sibHalf = ""
  if(any(niece)) {
    j = which(slot == "SIB" & halfSib)
    if(length(j)) {
      sibHalf = halfSlot[j]
    } else if(any(halfNiece)) {
      sibHalf = if(side[which(halfNiece)[1]] == "P") "PHM" else "MHF"
    }
  }
  halfSlots = unique.default(c(halfSlot[halfSib], sibHalf[nzchar(sibHalf)]))

  # Only create ancestors and intermediate relatives actually needed
  anc = any(type %in% c("father", "mother", "parent", "brother", "sister", "sibling",
                        "grandfather", "grandmother", "uncle", "aunt", "niece", "nephew"))
  pat = any(slot %in% c("FF", "FM")) || any(uncle & side == "P")
  mat = any(slot %in% c("MF", "MM")) || any(uncle & side == "M")
  roles = c("MP", "F", "M", "FF", "FM", "MF", "MM", spSlots, "SIB", "SIBSP", "CH", "CHSP", halfSlots)
  active = c(TRUE, anc, anc, pat, pat, mat, mat, rep(TRUE, length(spSlots)), any(niece), any(niece),
             any(grandchild), any(grandchild), rep(TRUE, length(halfSlots)))
  names(active) = roles
  id = setNames(paste0(family, "_", roles), roles)
  id["MP"] = mp
  id[slot[nzchar(slot)]] = a$label[nzchar(slot)]

  sibSex = if(any(slot == "SIB") && type[slot == "SIB"] == "sister") 2L else 1L
  chSex = if(any(slot == "CH") && type[slot == "CH"] == "daughter") 2L else 1L
  sex = c(MP = mpSex, F = 1L, M = 2L, FF = 1L, FM = 2L, MF = 1L, MM = 2L,
          setNames(rep(if(mpSex == 2) 1L else 2L, length(spSlots)), spSlots),
          SIB = sibSex, SIBSP = 3L - sibSex,
          CH = chSex, CHSP = 3L - chSex,
          setNames(ifelse(startsWith(halfSlots, "PHM"), 2L, 1L), halfSlots))

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
    fa["SIB"] = if(startsWith(sibHalf, "MHF")) sibHalf else "F"
    mo["SIB"] = if(startsWith(sibHalf, "PHM")) sibHalf else "M"
  }
  if(any(grandchild)) {
    chSP = if(any(slot == "CH")) spouseSlot[slot == "CH"] else spSlots[1]
    fa["CH"] = if(mpSex == 1) "MP" else chSP
    mo["CH"] = if(mpSex == 1) chSP else "MP"
  }

  # Parent links for references not occupying a structural position
  rf = rm = rep("", n)
  rf[sibling] = "F"
  rm[sibling] = "M"
  rf[halfSib & side == "P"] = "F"
  rm[halfSib & side == "P"] = halfSlot[halfSib & side == "P"]
  rf[halfSib & side == "M"] = halfSlot[halfSib & side == "M"]
  rm[halfSib & side == "M"] = "M"
  rf[uncle & side == "P"] = "FF"
  rm[uncle & side == "P"] = "FM"
  rf[uncle & side == "M"] = "MF"
  rm[uncle & side == "M"] = "MM"
  rf[child] = if(mpSex == 1) "MP" else spouseSlot[child]
  rm[child] = if(mpSex == 1) spouseSlot[child] else "MP"
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
