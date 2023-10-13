#' Nonidentifiable missing persons
#'
#' A missing person in a DVI case is *nonidentifiable* if unrelated to all
#' (genotyped) reference individuals and all other missing persons in the
#' reference family. It is often wise to ignore such individuals in [jointDVI()]
#' and other analyses, to relieve the computational burden.
#'
#' The implementation uses `ribd::kinship()` to identify individuals having
#' kinship coefficient 0 with all relevant individuals.
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
  
  k = ribd::kinship(dvi$am)
  
  # Individuals to consider: references and missing
  typed = typedMembers(dvi$am)
  ids = c(typed, dvi$missing)
  
  # Subset of kinship matrix
  kk = k[ids, ids, drop = FALSE]
  
  # Hack to simplify next step: Set diagonal entries to 0
  diag(kk) = 0
  
  # Columns with only zeroes = nonidentifiable
  nonident = colnames(kk)[colSums(kk) == 0]
  
  if(!length(nonident))
    return(list(nonidentifiable = nonident,
                dviReduced = dvi,
                summary = NULL))
  
  # Remove from `missing` slot
  if(length(nonident))
    dviRed = subsetDVI(dvi, missing = setdiff(dvi$missing, nonident), verbose = FALSE)
  else
    dviRed = dvi
  
  # For summary: Family of nonidentifiable indivs
  fams = getFamily(dvi, ids = nonident)
  
  # For summary: String of refs/missing for each
  cmts = lapply(seq_along(nonident), function(i) {
    idsi = intersect(ids, labels(dvi$am[[i]])) |> setdiff(nonident[i])
    toString(idsi)
  })
  
  # Build summary report
  summary = data.frame(Family = fams,
                       Missing = nonident,
                       Conclusion = "nonidentifiable",
                       Comment = paste("unrelated to", unlist(cmts, use.names = FALSE)),
                       row.names = NULL)
                        
  list(nonidentifiable = nonident,
       dviReduced = dviRed,
       summary = summary)
}
