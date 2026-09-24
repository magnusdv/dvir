#' Identify and merge matching PM samples
#'
#' Computes the direct-matching LR of each pair of samples, and merges the matching
#' samples.
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
#' @param method A keyword indicating how to merge matching samples. See Details.
#' @param dropout Allelic dropout rate. Default: 0.
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
                   dropout = 0, verbose = TRUE) {
  
  n = length(pm)
  method = match.arg(method)
  
  if(!is.numeric(threshold) || length(threshold) != 1 || is.na(threshold) || threshold <= 0)
      stop2("`threshold` must be a positive number")

  if(!is.numeric(dropout) || length(dropout) != 1 || is.na(dropout) || dropout < 0 || dropout >= 1)
      stop2("`dropout` must be a number in [0, 1)")

  if(verbose) {
    msg = c(sprintf("Number of singletons: %d", n),
            sprintf("LR threshold: %g", threshold),
            sprintf("Allelic dropout rate: %g", dropout),
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
  
  # Precompute likelihood of each singleton
  liks = lapply(pm, likelihood, dropout = dropout)
  
  # LR matrix (upper triangular)
  LRs = matrix(0, nrow = n, ncol = n, dimnames = list(ids, ids))
  for(i in 1:(n-1)) for(j in (i+1):n)
    LRs[i,j] = directMatch(pm[[i]], pm[[j]], 
                           g1 = g[i, ], g2 = g[j, ], 
                           dropout = dropout, 
                           .lik1 = liks[[i]], .lik2 = liks[[j]],
                           .skipChecks = TRUE)
  
  # Find connected groups of matching samples
  clust = list()
  for(i in 1:n) {
    rmatch = c(i, which(LRs[i, ] >= threshold))
    hit = which(vapply(clust, function(z) any(rmatch %in% z), logical(1)))

    if(!length(hit)) # create new comp
      clust[[length(clust) + 1]] = rmatch 
    else {
      clust[[hit[1]]] = unique.default(c(rmatch, unlist(clust[hit], use.names = FALSE)))
      clust[hit[-1]] = NULL
    }
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
#' Computes the likelihood ratio comparing if two samples are from the same individual or
#' from unrelated individuals.
#'
#' @param x,y Typed singletons.
#' @param g1,g2 (Optional) Named character vectors with genotypes for `x` and `y`
#'   respectively.
#' @param dropout Allelic dropout rate. Default: 0.
#' @param .lik1,.lik2 (For internal use.) Precomputed likelihoods for `x` and `y`
#'   respectively.
#' @param .skipChecks A logical indicating that various input checks can be skipped, e.g.
#'   when called by `mergePM()`.
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
directMatch = function(x, y, g1 = NULL, g2 = NULL, dropout = 0, 
                       .lik1 = NULL, .lik2 = NULL, .skipChecks = FALSE) {
  if(!.skipChecks) {
    if(!is.singleton(x))
      stop2("First argument is not a singleton: ", class(x)[1])
    if(!is.singleton(y))
      stop2("Second argument is not a singleton: ", class(y)[1])
    
    if(!is.numeric(dropout) || length(dropout) != 1 || is.na(dropout) || dropout < 0 || dropout >= 1)
      stop2("`dropout` must be a number in [0, 1)")

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
    
    x = selectMarkers(x, commonM)
    y = selectMarkers(y, commonM)
  
    g1 = g1[commonM]
    g2 = g2[commonM]
  }
  
  if(x$SEX != y$SEX && x$SEX * y$SEX > 0)
    return(0)
  
  miss1 = g1 == "-/-"
  miss2 = g2 == "-/-"
  
  nonmiss = which(!miss1 & !miss2)
  if(!length(nonmiss))
    return(1)
  
  if(dropout == 0) {
    if(!all(miss1 | miss2 | g1 == g2))
      return(0)
  
    lik = .lik1[nonmiss] %||% likelihood(x, markers = nonmiss)
    return(prod(1/lik))
  }
  
  # With dropout
  lik1 = .lik1[nonmiss] %||% likelihood(x, markers = nonmiss, dropout = dropout)
  lik2 = .lik2[nonmiss] %||% likelihood(y, markers = nonmiss, dropout = dropout)
  
  numer = vapply(nonmiss, \(i)
    .jointMatchLik(x$MARKERS[[i]], y$MARKERS[[i]], dropout), numeric(1))

  prod(numer/(lik1 * lik2))
}


.jointMatchLik = function(m1, m2, dropout) {
  a = m1[1, 1]; b = m1[1, 2]
  c = m2[1, 1]; d = m2[1, 2]
  p = attr(m1, "afreq")
  s = 1 - dropout

  het1 = a != b
  het2 = c != d

  if(het1) {
    hw = 2 * p[a] * p[b]
    if(het2)
      return(hw * s^4 * ((a == c && b == d) || (a == d && b == c)))
    return(hw * dropout * s^3 * (c == a || c == b))
  }

  if(het2)
    return(2 * p[c] * p[d] * dropout * s^3 * (a == c || a == d))

  if(a != c)
    return(2 * p[a] * p[c] * (dropout * s)^2)

  pa = p[a]
  pa^2 * (1 - dropout^2)^2 + 2 * pa * (1 - pa) * (dropout * s)^2
}