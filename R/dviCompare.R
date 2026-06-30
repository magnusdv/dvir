#' Compare DVI approaches
#'
#' Legacy helper for comparing different computational approaches to DVI.
#'
#' NOTE: This is a legacy function, mainly kept for reproducibility of older
#' analyses. It is not part of the recommended `dvir` workflow. Use with caution.
#'
#' The following methods are available through the `methods` parameter:
#'
#' 1. Sequential, without LR updates.
#'
#' 2. Sequential, with LR updates.
#'
#' 3. Sequential (undisputed) + joint (remaining). Always returns the most
#' likely solution(s).
#'
#' 4. Joint - brute force. Always returns the most likely solution(s).
#'
#' 5. Like 3, but returns the winner(s) only if LR > `threshold`; otherwise the
#' empty assignment.
#'
#' 6. Like 4, but returns the winner(s) only if LR > `threshold`; otherwise the
#' empty assignment.
#'
#' @param dvi A `dviData` object, typically created with [dviData()]. If
#'   `simulate = FALSE`, a list of such objects.
#' @param true A character of the same length as `dvi$pm`, with the true
#'   solution, e.g. `true = c("M2", "M3", "*")` if the truth is V1 = M2, V2 =
#'   M3 and V3 unmatched.
#' @param refs Character vector with names of the reference individuals. By
#'   default the typed members of `dvi$am`.
#' @param methods A subset of the numbers 1,2,3,4,5,6.
#' @param markers If `simulate = FALSE`: A vector indicating which markers
#'   should be used.
#' @param threshold An LR threshold passed on to the sequential methods.
#' @param simulate A logical, indicating if simulations should be performed.
#' @param db A frequency database used for simulation, e.g.
#'   `forrel::NorwegianFrequencies`. By default the frequencies attached to
#'   `dvi$am` are used.
#' @param Nsim A positive integer; the number of simulations.
#' @param returnSims A logical: If TRUE, the simulated data are returned without
#'   any DVI comparison.
#' @param seed A seed for the random number generator, or NULL.
#' @param numCores Deprecated and ignored.
#' @param verbose A logical.
#'
#' @return A list of solution frequencies for each method, and a vector of true
#'   positive rates for each method.
#'
#' @examples
#' db = forrel::NorwegianFrequencies[1:6]
#'
#' # True solution: V1 = M1, while V2/V3 are the siblings M2/M3
#' true = c("M1", "M2", "M3")
#'
#' \donttest{
#' # Increase Nsim for accurate estimates
#' dviCompare(fire, true = true, db = db, Nsim = 2, seed = 1)
#'
#' # Simulations can also be generated first
#' sims = dviCompare(fire, true = true, db = db, Nsim = 2, seed = 1,
#'                   returnSims = TRUE)
#' dviCompare(sims, true = true, simulate = FALSE)
#' }
#' 
#' @importFrom forrel profileSim
#' @export
dviCompare = function(dvi, true, refs = NULL, methods = 1:6, 
                      markers = NULL, threshold = 1, simulate = TRUE, 
                      db = NULL, Nsim = 1, returnSims = FALSE, 
                      seed = NULL, numCores = 1, verbose = TRUE) {
  
  if(verbose)
    message("NOTE: `dviCompare()` is a legacy function. Use with caution.")
  
  if(!identical(numCores, 1))
    warning("`numCores` is deprecated and currently ignored", call. = FALSE)
  
  st = Sys.time()
  threshold = max(threshold, 1 + .Machine$double.eps)
  
  dvi1 = if(simulate) dvi else if(inherits(dvi, "dviData")) dvi else dvi[[1]]
  dvi1 = consolidateDVI(dvi1)
  
  refs = if(is.null(refs)) typedMembers(dvi1$am) else as.character(refs)
  true = as.character(true)
  nsim = if(simulate) Nsim else if(inherits(dvi, "dviData")) 1 else length(dvi)
  
  if(verbose) {
    print(dvi1)
    cat("\nParameters for DVI comparison:\n")
    cat(" True solution:", toString(true), "\n")
    cat(" Simulate data:", simulate, "\n")
    cat(" Number of sims:", nsim, "\n")
    cat(" Reference IDs:", toString(refs), "\n")
    cat(" LR threshold:", threshold, "\n\n")
  }
  
  if(simulate) {
    dvi = dvi1
    
    pm = dvi$pm
    am = dvi$am
    missing = as.character(dvi$missing)
    
    vics = names(pm)
    isMatch = true != "*"
    
    stopifnot(length(true) == length(vics), 
              all(true[isMatch] %in% missing), 
              all(getSex(am, true[isMatch]) == getSex(pm)[isMatch]))
    
    if(!is.null(seed))
      set.seed(seed)
    
    db = db %||% getDatabase(dvi, check = FALSE)
    am = setMarkers(am, locusAttributes = db)
    pm = setMarkers(pm, locusAttributes = db)
    
    # Simulate AM
    AMsims = profileSim(am, N = Nsim, ids = c(refs, true[isMatch]), simplify1 = FALSE)
    
    # Simulate the unrelated victims
    PMsims = profileSim(pm, N = Nsim, ids = vics[!isMatch], simplify1 = FALSE)
    
    # For the true matches, transfer from MP to vics
    PMsims = lapply(seq_len(Nsim), function(i) 
      transferMarkers(from = AMsims[[i]], to = PMsims[[i]], 
                      idsFrom = true[isMatch], idsTo = vics[isMatch], erase = FALSE))
    
    # Remove data from missing
    AMsims = lapply(AMsims, function(s) removeGenotypes(s, missing))
    
    # Collect into list of dviData objects
    dviSims = lapply(seq_len(Nsim), function(i) 
      dviData(pm = PMsims[[i]], am = AMsims[[i]], missing = missing))
    
    # Return sims
    if(returnSims)
      return(dviSims)
  }
  else {
    dviSims = if(inherits(dvi, "dviData")) list(dvi) else dvi
    dviSims = lapply(dviSims, consolidateDVI)
    
    vics = names(dvi1$pm)
    stopifnot(length(true) == length(vics), 
              all(true %in% c(dvi1$missing, "*")),
              .mysetequal(refs, typedMembers(dvi1$am)))
    
    if(!is.null(markers))
      stop2("Marker selection not implemented.")
  }
  
  N = length(dviSims)
  
  # DVI functions
  funs = list(
    function(x) sequentialDVI(x, threshold = threshold, updateLR = FALSE, check = FALSE, verbose = FALSE)$matches,
    function(x) sequentialDVI(x, threshold = threshold, updateLR = TRUE, check = FALSE, verbose = FALSE)$matches,
    function(x) top(jointDVI(x, undisputed = TRUE, threshold = threshold, check = FALSE, verbose = FALSE)),
    function(x) top(jointDVI(x, undisputed = FALSE, check = FALSE, verbose = FALSE)),
    function(x) top(jointDVI(x, undisputed = TRUE, threshold = threshold, check = FALSE, verbose = FALSE), threshold),
    function(x) top(jointDVI(x, undisputed = FALSE, check = FALSE, verbose = FALSE), threshold)
  )
  
  res = lapply(methods, function(m) {
    z = summar(lapply(dviSims, funs[[m]]))
    if(verbose)
      print(setNames(list(z), paste0("method", m)))
    z
  })
  
  names(res) = paste0("method", methods)
  
  # True positive rates
  TPR = sapply(res, function(r) as.numeric(r[paste(true, collapse = "-")]))
  TPR[is.na(TPR)] = 0
  res$TPR = TPR
  
  if(verbose)
    cat("Total time used:", format(Sys.time() - st, digits = 3), "\n")
  
  res
}

# Utility for summarising list of solutions
#' @importFrom stats aggregate
summar = function(x) {
  sols = lapply(x, apply, 1, paste, collapse = "-")
  L = lengths(sols)
  df = data.frame(sol = unlist(sols), wei = rep(1/L, L))
  agg = aggregate(wei ~ sol, data = df, FUN = sum)
  freqs = agg$wei/length(x)
  names(freqs) = agg$sol
  sort(freqs, decreasing = TRUE)
}

# Utility for including ties in output of jointDVI()
top = function(res, thresh = 1) {
  nvic = match("loglik", colnames(res), nomatch = 1L) - 1
  if(res$LR[1] < thresh) return(rbind(rep("*", nvic)))
  mx = res$LR == res$LR[1]
  res[mx, seq_len(nvic), drop = FALSE]
}