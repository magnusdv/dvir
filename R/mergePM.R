#' Identify and merge matching PM samples
#'
#' Computes pairwise direct-match LRs for post-mortem samples, identifies groups of
#' matching samples, and reduces each group to a single profile.
#'
#' Groups are defined as connected components of pairs whose LR is at least `threshold`.
#' Thus, samples may belong to the same group even if their pairwise LR is below the
#' threshold, provided they are connected through other samples.
#'
#' The available merging methods are:
#'
#' * `"combine"`: Combine the observed alleles across all samples in the group.
#' With dropout, discrepancies compatible with allelic dropout are allowed. Without
#' dropout, complete genotypes must agree. Markers that cannot be combined are set to
#' missing and reported in `problems`.
#' 
#' * `"first"`: Retain the first sample in each group, according to input order.
#'
#' * `"mostcomplete"`: Retain the sample with the most non-missing genotypes.
#' 
#' The names of the resulting clusters are controlled by `names`: `"first"` uses
#' the first sample, `"mostcomplete"` the most complete sample, while `"combine"` joins
#' all sample names with `"_"`.
#'
#' @param pm A list of typed singletons.
#' @param threshold LR threshold for positive identification.
#' @param method A keyword indicating how matching samples should be merged. See Details.
#' @param names A keyword controlling the names of merged samples; one of
#'   `"combine"`, `"first"` or `"mostcomplete"`.
#' @param dropout Allelic dropout probability. Default: 0.
#' @param verbose A logical.
#'
#' @returns A list with the following entries:
#'
#' * `groups`: The groups of matching samples.
#'
#' * `LRmat`: A symmetric matrix containing all pairwise direct-match LRs.
#'
#' * `nonmissing`: The number of non-missing genotypes for each sample.
#'
#' * `pmReduced`: The reduced list of PM samples.
#'
#' * `problems`: For `method = "combine"`, a named list of markers that could
#'   not be combined and were set to missing. Empty otherwise.
#'
#' @references Dørum G, Kling D, Baeza-Richer C, García-Magariños M, Sæbø S, Desmyter S,
#'   Egeland T (2015). "Models and implementation for relationship problems with dropout".
#'   *International Journal of Legal Medicine*, 129, 411-423.
#'   \doi{10.1007/s00414-014-1046-5}
#'
#' @seealso [directMatch()].
#'
#' @examples
#' afr = c("1" = 0.1, "2" = 0.9)
#'
#' pm = singletons(c("V1", "V2", "V3")) |>
#'   addMarker(V1 = "1/1", V2 = "1/1", V3 = "2/2",
#'             afreq = afr, name = "M1") |>
#'   addMarker(V1 = NA, V2 = "2/2", V3 = "1/2",
#'             afreq = afr, name = "M2")
#'
#' mergePM(pm, threshold = 10, verbose = FALSE)
#' mergePM(pm, threshold = 10, method = "mostcomplete", verbose = FALSE)
#'
#' @export
mergePM = function(pm, threshold = 1e4,
                   method = c("combine", "first", "mostcomplete"), 
                   names = c("combine", "first", "mostcomplete"), 
                   dropout = 0, verbose = TRUE) {
  
  if(!all(vapply(pm, is.singleton, logical(1))))
    stop2("First argument must be a list of singletons")
  
  n = length(pm)
  method = match.arg(method)
  names = match.arg(names)
  
  if(!is.numeric(threshold) || length(threshold) != 1 || is.na(threshold) || threshold <= 0)
      stop2("`threshold` must be a positive number")

  if(!is.numeric(dropout) || length(dropout) != 1 || is.na(dropout) || dropout < 0 || dropout >= 1)
      stop2("`dropout` must be a number in [0, 1)")

  if(verbose) {
    msg = c(sprintf("Number of singletons: %d", n),
            sprintf("LR threshold: %g", threshold),
            sprintf("Dropout prob: %g", dropout),
            sprintf("Merging method: %s", method),
            sprintf("Naming method: %s", names))
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

  # Generate names for the cluster groups
  gnames = switch(names,
    first = vapply(groups, \(g) g[1], character(1)),
    mostcomplete = vapply(groups, \(g) g[which.max(nonmissing[g])], character(1)),
    combine = vapply(groups, \(g) paste(g, collapse = "_"), character(1))
  )
  
  # Put most complete sample first when relevant
  if(method %in% c("mostcomplete", "combine"))
    groups = lapply(groups, function(g) g[order(nonmissing[g], decreasing = TRUE)])
  
  names(groups) = gnames
  
  # Merge matching samples
  problems = list()
  
  if(method == "combine") {
    comb = lapply(groups, \(g) .combinePM(pm[g], withDropout = dropout > 0))
  
    pmReduced = lapply(comb, `[[`, "profile")
    problems = lapply(comb, `[[`, "problems")
    problems = problems[lengths(problems) > 0]
  }
  else {
    # Keep first (also works for mostcomplete after sorting!)
    keep = vapply(groups, function(g) g[1], character(1))
    pmReduced = pm[keep]
  }
  
  names(pmReduced) = names(groups)
  
  # If needed, rename singletons internally also
  newlabs = names(pmReduced)
  if(!identical(labels(pmReduced), newlabs))
    pmReduced = relabel(pmReduced, new = newlabs)
  
  # Make LR matrix symmetric
  LRmat = LRs + t.default(LRs)
  
  if(verbose) {
    cat("-----\n")
    clust = groups[lengths(groups) > 1]
    ncl = length(clust)
    if(ncl) {
      ss = unlist(lapply(names(clust), function(nm) {
        prob = problems[[nm]]
        note = if(length(prob)) sprintf(" (inconsistent marker: %s)", toString(prob)) else ""
        sprintf(" * [%s]%s\n", toString(clust[[nm]]), note)
      }))
      cat(sprintf("%d cluster%s identified:\n%s", ncl, if(ncl != 1) "s" else "", ss), sep = "")
    }
    else
      cat("No clusters identified\n")
  }
  
  list(groups = groups,
       LRmat = LRmat,
       nonmissing = nonmissing,
       pmReduced = pmReduced,
       problems = problems)
}



#' Direct match LR
#'
#' Computes the likelihood ratio comparing the hypotheses that two samples originate from
#' the same individual or from two unrelated individuals.
#'
#' For a single marker, the LR is computed as
#'
#' \deqn{LR = \frac{P(G_1,G_2 \mid H_1)}
#'                  {P(G_1 \mid H_2)P(G_2 \mid H_2)},}
#'
#' where `G1` and `G2` are the observed genotypes of the two samples. `H1` states that
#' the samples originate from the same individual, and `H2` that they originate from
#' unrelated individuals. The overall LR is obtained by multiplying the marker-wise LRs.
#' 
#' With `dropout = 0`, discordant non-missing genotypes give LR = 0. For positive dropout
#' we use the model of Dørum et al. (2015), where alleles drop out independently with
#' probability `d`. In particular,
#'
#' \deqn{P(a/b \mid a/b) = (1-d)^2,} 
#' \deqn{P(a/a \mid a/b) = d(1-d),} 
#' \deqn{P(a/a \mid a/a) = 1-d^2.}
#'
#' Thus, apparently discordant genotypes may have a positive LR when explained by allelic
#' dropout. A marker missing in either sample contributes LR = 1.
#'
#' @param x,y Typed singletons.
#' @param g1,g2 Optional named character vectors containing precomputed genotypes for `x`
#'   and `y`, respectively.
#' @param dropout Allelic dropout probability. Default: 0.
#' @param .lik1,.lik2 For internal use; precomputed likelihoods for `x` and `y`.
#' @param .skipChecks For internal use; skip input checks.
#'
#' @return A single number.
#'
#' @references Dørum et. al (2015). "Models and implementation for relationship problems
#'   with dropout". *International Journal of Legal Medicine*, 129, 411-423.
#'   \doi{10.1007/s00414-014-1046-5}
#'
#' @seealso [mergePM()].
#'
#' @examples
#' afr = c("1" = 0.1, "2" = 0.9)
#' pm = singletons(c("V1", "V2", "V3")) |>
#'   addMarker(V1 = "1/2", V2 = "1/2", V3 = "1/1",
#'             afreq = afr, name = "M")
#'
#' directMatch(pm[[1]], pm[[2]])
#' directMatch(pm[[1]], pm[[3]])
#' directMatch(pm[[1]], pm[[3]], dropout = 0.1)
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


.combinePM = function(pm, withDropout = FALSE) {
  if(length(pm) == 1)
    return(list(profile = pm[[1]], problems = character()))

  z = pm[[1]]
  nM = length(z$MARKERS)

  # Allele codes: 2 x markers x samples
  a = vapply(pm, function(x) unlist(x$MARKERS, use.names = FALSE), numeric(2 * nM))
  dim(a) = c(2, nM, length(pm))

  problems = character()

  for(m in seq_len(nM)) {
    am = a[, m, ]
    u = unique.default(am[am > 0])
    bad = length(u) > 2

    # Without dropout, complete genotypes must agree
    if(!bad && !withDropout) {
      full = colSums(am > 0) == 2
      if(any(full)) {
        g = am[, full, drop = FALSE]
        g0 = g[, 1]
        same = (g[1, ] == g0[1] & g[2, ] == g0[2]) |
               (g[1, ] == g0[2] & g[2, ] == g0[1])
        bad = !all(same) || anyNA(match(u, g0))
      }
    }

    if(bad) {
      z$MARKERS[[m]][1, ] = 0L
      problems = c(problems, attr(z$MARKERS[[m]], "name"))
      next
    }

    # Two observed alleles imply heterozygosity
    if(length(u) == 2)
      geno = u

    # One observed allele: retain a/a if seen, otherwise a/-
    else if(length(u) == 1) {
      hom = any(am[1, ] == u & am[2, ] == u)
      geno = if(hom) rep.int(u, 2) else c(u, 0L)
    }

    # No observed alleles
    else
      geno = c(0L, 0L)

    z$MARKERS[[m]][1, ] = geno
  }

  list(profile = z, problems = problems)
}