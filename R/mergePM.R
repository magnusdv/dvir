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
#' @param verbose A logical.

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
#' pm = singletons(c("V1", "V2", "V3")) |> 
#'   addMarker(V1 = "1/1", V2 = "2/2", V3 = "1/1", 
#'             afreq = c("1" = 0.01, "2" = 0.99), name = "L1")
#' 
#' mergePM(pm)
#'
#' @export
mergePM = function(pm, threshold = 1e4, method = c("mostcomplete", "first", "combine"), 
                   verbose = TRUE) {
  
  n = length(pm)
  method = match.arg(method)
  
  if(verbose) {
    msg = c(sprintf("Number of singletons: %d", n),
            sprintf("LR threshold: %g", threshold),
            sprintf("Merging method: '%s'", method))
    cat(msg, sep = "\n")
  }
  if(n < 2) {
    if(verbose) cat("Nothing to do\n")
    return(list())
  }
  
  g = getGenotypes(pm)
  ids = rownames(g)
  names(pm) = ids
  
  # Number of non-missing for each sample
  nonmissing = rowSums(g != "-/-")
  
  # LR matrix (upper triangular)
  LRs = matrix(0, nrow = n, ncol = n, dimnames = list(ids, ids))
  for(i in 1:(n-1)) for(j in (i+1):n)
    LRs[i,j] = directMatch(pm[[i]], pm[[j]], g1 = g[i, ], g2 = g[j, ], 
                           .skipChecks = TRUE)
  
  # Find clusters of matching samples. NB: Indices!
  clust = list()
  for(i in 1:n) {
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
  
  if(verbose) {
    cat("-----\n")
    clust = groups[lengths(groups) > 1]
    if(length(clust)) {
      s = unlist(lapply(clust, function(g) sprintf(" * [%s]\n", toString(g))), use.names = FALSE)
      cat("Groups of matching samples:\n", s, sep = "")
    }
    else
      cat("Groups of matching samples: None\n")
  }
  
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
#' @param g1,g2 (Optional) Named character vectors with genotypes for `x` and
#'   `y` respectively.
#' @param .skipChecks A logical indicating that various input checks can be skipped,
#'   e.g. when called by `mergePM()`.
#'
#' @return A nonnegative likelihood ratio.
#' @seealso [mergePM()].
#'
#' @examples
#'
#' pm = singletons(c("V1", "V2", "V3")) |>
#'   addMarker(V1 = "1/1", V2 = "2/2", V3 = "1/1",
#'             afreq = c("1" = 0.01, "2" = 0.99), name = "L1")
#'
#' directMatch(pm[[1]], pm[[2]])
#' directMatch(pm[[1]], pm[[3]])
#'
#' @export
directMatch = function(x, y, g1 = NULL, g2 = NULL, .skipChecks = FALSE) {
  if(!.skipChecks) {
    if(!is.singleton(x))
      stop2("First argument is not a singleton: ", class(x)[1])
    if(!is.singleton(y))
      stop2("Second argument is not a singleton: ", class(y)[1])
    
    if(is.null(g1))
      g1 = getGenotypes(x)[1,]
    if(is.null(g2))
      g2 = getGenotypes(y)[1,]
    
    # Add missing name (occurs in cases with only 1 marker)
    if(is.null(names(g1)))
      names(g1) = name(x)
    if(is.null(names(g2)))
      names(g2) = name(y)
    
    commonM = .myintersect(names(g1), names(g2))
    if(!length(commonM)) {
      message("No shared markers")
      return(1)
    }
    
    g1 = g1[commonM]
    g2 = g2[commonM]
  }
  
  if(x$SEX != y$SEX && x$SEX * y$SEX > 0)
    return(0)
  
  miss1 = g1 == "-/-"
  miss2 = g2 == "-/-"
  if(!all(miss1 | miss2 | g1 == g2))
    return(0)
  
  nonmiss = which(!miss1 & !miss2)
  
  lik = likelihood(x, markers = nonmiss)
  prod(1/lik)
}


