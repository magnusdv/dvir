#' Identity and merge matching PM samples
#'
#' Computes the direct matching LR of each pair of samples, and merges the
#' matching samples.
#'
#' The available methods for merging matched samples are:
#'
#' * "mostcomplete": Use the sample with the highest number of non-missing genotypes
#'
#' * "first": Use the first in each group, according to the input order
#'
#' * "combine": Not implemented yet.
#'
#' @param pm A list of typed singletons.
#' @param threshold LR threshold for positive identification.
#' @param method A keyword indicating how to merging matching samples. See Details.
#'
#' @seealso [directMatch()].
#'
#' @returns A list with the following entries:
#' * `groups`: A list containing the groups of matching samples.
#'
#' * `LRmat`: A symmetric matrix (with 0s on the diagonal) containing the direct
#'   matching LR values.
#'
#' * `nonmissing`: A named vector reporting the number of non-missing genotypes 
#' for each sample.
#' 
#' * `pmReduced`: A list of singletons. If `use` is "best" or "first", this is 
#' a subset of the input `pm`.
#' 
#' @examples
#'
#' library(forrel)
#'
#' x = singleton("a") |> profileSim(markers = NorwegianFrequencies[1:4])
#' y = relabel(x, new = "b")
#'
#' # Remove one genotype to make the example more interesting
#' x = setGenotype(x, id = "a", marker = 1, geno = "-/-")
#'
#' mergePM(list(x, y))
#'
#' @export
mergePM = function(pm, threshold = 1e4, method = c("mostcomplete", "first", "combine")) {
  n = length(pm)
  if(n < 2)
    return(list())
  
  method = match.arg(method)
  g = getGenotypes(pm)
  ids = rownames(g)
  names(pm) = ids
  
  # Number of non-missing for each sample
  nonmissing = rowSums(g != "-/-")
  
  # LR matrix (upper triangular)
  LRs = matrix(0, nrow = n, ncol = n, dimnames = list(ids, ids))
  for(i in 1:(n-1)) for(j in (i+1):n)
    LRs[i,j] = directMatch(pm[[i]], pm[[j]], geno1 = g[i, ], geno2 = g[j, ])
  
  # Find clusters of matching samples. NB: Indices!
  clust = list()
  for(i in 1:(n-1)) {
    
    # Matches in row i (including diagonal entry)
    rmatch = c(i, which(LRs[i, ] >= threshold))
    new = TRUE
    
    # Loop over current clusters
    K = length(clust)
    for(k in seq_len(K)) {
      this = clust[[k]]
      if(any(rmatch %in% this)) {
        clust[[k]] = c(this, rmatch)
        new = FALSE
        break
      }
    }
    # If not in previous comps, create new
    if(new)
      clust[[K+1]] = rmatch
  }
  
  # Convert indices to names (sorted by input order)
  groups = lapply(clust, function(idx) ids[sort.default(unique.default(idx))])
  
  # For "mostcomplete", re-sort and add names
  if(method == "mostcomplete") {
    groups = lapply(groups, function(g) g[order(nonmissing[g], decreasing = TRUE)])
    names(groups) = sapply(groups, '[', 1)
  }
    
  # Merge matching samples
  pmReduced = switch(method,
    mostcomplete = pm[names(groups)],
    first = pm[unlist(lapply(groups, function(g) g[1]))],
    combine = stop2("Method 'combine' is not implemented yet")
  )
    
  # Make LR matrix symmetric (with 0 on diag)
  LRmat = LRs + t.default(LRs)
  
  list(groups = groups, 
       LRmat = LRmat,
       nonmissing = nonmissing,
       pmReduced = pmReduced)
}



#' Direct match LR
#'
#' Computes the likelihood ratio comparing if two samples are from the same
#' individual or from unrelated individuals.
#'
#' @param x,y Typed singletons.
#' @param geno1,geno2 (Optional) Named character vectors with genotypes for `x`
#'   and `y` respectively.
#'
#' @return A nonnegative likelihood ratio.
#' @seealso [mergePM()].
#' 
#' @examples
#' library(forrel)
#' 
#' x = singleton("a") |> profileSim(markers = NorwegianFrequencies[1:4])
#' y = relabel(x, new = "b")
#'
#' # Remove one genotype to make the example more interesting
#' x = setGenotype(x, id = "a", marker = 1, geno = "-/-")
#'
#' directMatch(x, y)
#'
#' @export
directMatch = function(x, y, geno1 = NULL, geno2 = NULL) {
  if(!is.singleton(x))
    stop2("First argument is not a singleton: ", class(x)[1])
  if(!is.singleton(y))
    stop2("Second argument is not a singleton: ", class(y)[1])
  
  if(x$SEX != y$SEX && x$SEX * y$SEX > 0)
    return(0)
  
  g1 = geno1 %||% getGenotypes(x)[1,]
  g2 = geno2 %||% getGenotypes(y)[1,]
  
  commonM = intersect(names(g1), names(g2))
  if(!length(commonM)) {
    message("No shared markers")
    return(1)
  }
  
  g1 = g1[commonM]
  g2 = g2[commonM]
  miss1 = g1 == "-/-"
  miss2 = g2 == "-/-"
  if(!all(miss1 | miss2 | g1 == g2))
    return(0)
  
  nonmiss = commonM[!miss1 & !miss2]
  
  lik = likelihood(x, markers = nonmiss)
  prod(1/lik)
}


