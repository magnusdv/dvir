#' Nonidentifiable missing persons
#'
#' A missing person in a DVI case is nonidentifiable if they are not related to
#' any (genotyped) reference individuals, neither directly nor through other
#' missing persons in the family. This function identifies such individuals, and
#' reduces the DVI dataset accordingly.
#'
#' Two individuals are unrelated if their kinship coefficient is zero. The
#' implementation uses `ribd::kinship()` to calculate all pairwise kinship
#' coefficients.
#'
#' @param dvi A `dviData` object, typically created with [dviData()].
#'
#' @return A list with the following entries:
#'   * `nonidentifiable`: A character vector (possibly empty) with the names of
#'   the nonidentifiable missing persons.
#'   * `dviReduced`: A reduced `dviData` object, where the nonidentifiable
#'   individuals are removed from the list of missing persons. If there are no
#'   `nonidentifiable`, this is just a copy of `dvi`.
#'   * `summary`: A data frame summarising the findings.
#'
#' @examples
#' # Example 1: No nonidentifiables in dataset `example1`
#' findNonidentifiable(example1)
#'
#' # Example 2: Add nonidentifiable person "A"
#' amNew = example1$am[[1]] |>
#'   addSon(parents = c("NN", "A"))
#' missNew = c(example1$missing, "A")
#'
#' dvi = dviData(pm = example1$pm, am = amNew, missing = missNew)
#' plotDVI(dvi, textAbove = c(A = "nonidentif."))
#'
#' findNonidentifiable(dvi)
#'
#' @export
findNonidentifiable = function(dvi) {
  dvi = consolidateDVI(dvi)
  miss = dvi$missing
  refs = typedMembers(dvi$am)
  
  # Kinship matrix of refs and missing
  ids = c(typedMembers(dvi$am), dvi$missing)
  k = ribd::kinship(dvi$am, ids = ids, simplify = FALSE) # TODO: across = FALSE ?
  
  # Columns with only zeroes = nonidentifiable
  nonident = .nonidentif(k, miss, refs)
  
  # If none, return unchanged
  if(!length(nonident))
    return(list(nonidentifiable = nonident,
                dviReduced = dvi,
                summary = NULL))
  
  # Reduce by removing nonidentifiable from "missing"
  dviRed = subsetDVI(dvi, missing = setdiff(miss, nonident), verbose = FALSE)
  
  # For summary: Family of nonidentifiable indivs
  fams = getFamily(dvi, ids = nonident)

  # Special comment if there are no refs in a family
  nRefs = sapply(dvi$am, function(x) length(typedMembers(x)))

  # Comment column
  cmt = ifelse(nRefs[fams] == 0, "No references in family", "Unrelated to all references")
  
  # Build summary report
  summary = data.frame(Family = fams,
                       Missing = nonident,
                       Conclusion = "Nonidentifiable",
                       Comment = cmt,
                       row.names = NULL)
                        
  list(nonidentifiable = nonident,
       dviReduced = dviRed,
       summary = summary)
}

.nonidentif = function(k, miss, refs) {
  # Candidates for being nonidentifiable
  cand = miss
  
  # Individuals linked to the refs
  linked = refs
  
  # Iteratively remove those who are related to the previous links
  while(length(cand) > 0 && length(linked) > 0) {
    linked = names(which(rowSums(k[cand, linked, drop = FALSE]) > 0))
    cand = setdiff(cand, linked)
  }
  
  cand
}
