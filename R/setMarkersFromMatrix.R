#' Attach markers from an allele matrix
#'
#' A fast alternative to [pedtools::setMarkers()] for prepared allele matrices. All
#' pedigree components use the same alleles and frequencies, and untyped individuals
#' receive missing genotypes.
#'
#' @param x A `ped` object, a list of pedigrees, or a nested list of these.
#' @param alleleMatrix Matrix with two allele columns per marker and row names matching
#'   genotyped individuals. With `loci = NULL`, columns must be named in matching pairs,
#'   e.g. `D1, D1` or `D1 1, D1 2`.
#' @param loci Optional list of locus attributes, each containing `name`, `alleles`,
#'   `afreq`, and optionally `mutmod`. If NULL, observed alleles are assigned equal
#'   frequencies. All observed alleles must be included otherwise.
#'
#' @return The input pedigrees with markers attached, preserving list structure.
#'
#' @examples
#' peds = list(nuclearPed(), singleton("NN"))
#' db = forrel::NorwegianFrequencies[1:2]
#' ids = c("3", "NN")
#'
#' # Simulate profiles and extract alleles
#' x = forrel::profileSim(peds, ids = ids, markers = db, seed = 1)
#' a = getAlleles(x, ids = ids)
#' loc = getLocusAttributes(x)
#'
#' y1 = setMarkersFromMatrix(x, a)
#' y2 = setMarkersFromMatrix(x, a, loci = loc)
#'
#' stopifnot(identical(x, y2))
#'
#' # Note the difference:
#' afreq(y1, 1)
#' afreq(y2, 1)
#'
#' @export
setMarkersFromMatrix = function(x, alleleMatrix = NULL, loci = NULL) {
  if(is.null(alleleMatrix))
    return(x)
  if(is.null(rownames(alleleMatrix)))
    stop2("Allele matrix must have sample IDs as row names")

  # Derive shared locus attributes before processing individual pedigrees
  if(is.null(loci)) {
    nm = colnames(alleleMatrix)
    if(is.null(nm) || ncol(alleleMatrix) %% 2L)
      stop2("With loci = NULL, allele columns must be named in pairs")
    
    j = 2L*seq_len(ncol(alleleMatrix)/2L) - 1L
    same = nm[j] == nm[j + 1L]
    nm1 = ifelse(same, nm[j], sub("[-_. ]1$", "", nm[j]))
    nm2 = ifelse(same, nm[j + 1L], sub("[-_. ]2$", "", nm[j + 1L]))
    if(anyNA(nm1) || !identical(nm1, nm2) || any(!nzchar(nm1)))
      stop2("Allele columns must have matching marker names")
    
    loci = lapply(seq_along(j), function(k) {
      als = unique.default(as.character(alleleMatrix[, j[k] + 0:1]))
      als = sort(als[!is.na(als) & !als %in% c("", "-", "0", "NA")])
      if(!length(als)) 
        stop2("No observed alleles at marker: ", nm1[k])
      list(name = nm1[k], alleles = als, afreq = rep.int(1/length(als), length(als)))
    })
  }
  if(is.list(x) && !is.ped(x))
    return(lapply(x, setMarkersFromMatrix, alleleMatrix = alleleMatrix, loci = loci))
  if(!is.ped(x))
    stop2("Input must be a pedigree or a list of pedigrees")
  if(ncol(alleleMatrix) != 2L * length(loci))
    stop2("Number of allele columns must be twice the number of loci")

  # Match sample rows to pedigree members once
  ids = x$ID
  idx = match(rownames(alleleMatrix), ids, 0L)
  keep = idx > 0L
  idx = idx[keep]
  template = matrix(0L, nrow = length(ids), ncol = 2L)

  # Construct each marker using the supplied common allele set
  mlist = lapply(seq_along(loci), function(k) {
    loc = loci[[k]]
    m = template
    m[idx, ] = match(alleleMatrix[keep, (2L*k - 1L):(2L*k), drop = FALSE],
                      loc$alleles, nomatch = 0L)
    newMarker(m, alleles = loc$alleles, afreq = loc$afreq, name = loc$name,
              mutmod = loc$mutmod, pedmembers = ids, sex = x$SEX)
  })
  class(mlist) = "markerList"
  x$MARKERS = mlist
  x
}


ix = function(x, alleleMatrix = NULL, loci = NULL) {
  if(is.null(alleleMatrix)) 
    return(x)
  if(is.null(rownames(alleleMatrix)))
    stop2("Allele matrix must have sample IDs as row names")

  # Derive shared locus attributes before processing individual pedigrees
  if(is.null(loci)) {
    nm = colnames(alleleMatrix)
    if(is.null(nm) || ncol(alleleMatrix) %% 2L)
      stop2("With loci = NULL, allele columns must be named in pairs")
    
    j = 2L * seq_len(ncol(alleleMatrix)/2L) - 1L
    i = j[which(nm[j] != nm[j + 1L])]
    nm[i] = sub("[-_. ]1$", "", nm[i])
    nm[i + 1L] = sub("[-_. ]2$", "", nm[i + 1L])
    if(anyNA(nm) || !identical(nm[j], nm[j + 1L]) || any(!nzchar(nm[j])))
      stop2("Allele columns must have matching marker names")

    if(!is.character(alleleMatrix)) 
      storage.mode(alleleMatrix) = "character"
    
    miss = c("", "-", "0", "NA")
    loci = lapply(seq_along(j), function(k) {
      als = unique.default(alleleMatrix[, j[k] + 0:1])
      als = sort.int(als[!is.na(als) & !als %in% miss], method = "quick")
      nals = length(als)
      if(!nals)
        stop2("No observed alleles at marker: ", nm[j[k]])
      list(name = nm[j[k]], alleles = als, afreq = rep.int(1/nals, nals))
    })
  }
  if(ncol(alleleMatrix) != 2L * length(loci))
    stop2("Number of allele columns must be twice the number of loci")
  if(is.list(x) && !is.ped(x))
    return(lapply(x, setMarkersFromMatrix, alleleMatrix = alleleMatrix, loci = loci))
  if(!is.ped(x)) stop2("Input must be a pedigree or a list of pedigrees")

  # Match pedigree members and subset genotyped rows once
  ids = x$ID
  idx = match(ids, rownames(alleleMatrix))
  typed = which(!is.na(idx))
  a = alleleMatrix[idx[typed], , drop = FALSE]
  template = matrix(0L, nrow = length(ids), ncol = 2L)

  # Construct each marker using the supplied common allele set
  mlist = lapply(seq_along(loci), function(k) {
    loc = loci[[k]]
    m = template
    m[typed, ] = match(a[, (2L*k - 1L):(2L*k)], loc$alleles, nomatch = 0L)
    newMarker(m, alleles = loc$alleles, afreq = loc$afreq, name = loc$name,
              mutmod = loc$mutmod, pedmembers = ids, sex = x$SEX)
  })
  class(mlist) = "markerList"
  x$MARKERS = mlist
  x
}
