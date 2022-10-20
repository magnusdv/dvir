#' Compare DVI approaches
#'
#' Compare the efficiency of different computational approaches to DVI.
#'
#' The following methods are available for comparison, through the `methods`
#' parameter:
#'
#' 1. Sequential, without LR updates
#'
#' 2. Sequential, with LR updates
#'
#' 3. Sequential (undisputed) + joint (remaining). Always return the most likely
#' solution(s).
#'
#' 4. Joint - brute force. Always return the most likely solution(s).
#'
#' 5. Like 3, but return winner(s) only if LR > `threshold`; otherwise the empty
#' assignment.
#'
#' 6. Like 4, but return winner(s) only if LR > `threshold`; otherwise the empty
#' assignment.
#'
#' @param pm PM data: List of singletons
#' @param am AM data: A ped object or list of such.
#' @param missing Character vector with names of the missing persons.
#' @param true A character of the same length as `pm`, with the true solution,
#'   e.g., `true = c("M2", "*", "M3)` if the truth is V1 = M2 and V3 = M3.
#' @param refs Character vector with names of the reference individuals. By
#'   default the typed members of `am`.
#' @param methods A subset of the numbers 1,2,3,4,5,6.
#' @param markers If `simulate = FALSE`: A vector indicating which markers
#'   should be used.
#' @param threshold An LR threshold passed on to the sequential methods.
#' @param simulate A logical, indicating if simulations should be performed.
#' @param db A frequency database used for simulation, e.g.,
#'   `forrel::NorwegianFrequencies`. By default the frequencies attached to `am`
#'   are used.
#' @param Nsim A positive integer; the number of simulations.
#' @param returnSims A logical: If TRUE, the simulated data are returned without
#'   any DVI comparison.
#' @param seed A seed for the random number generator, or NULL.
#' @param numCores The number of cores used in parallelisation. Default: 1.
#' @param verbose A logical.
#'
#' @return A list of solution frequencies for each method, and a vector of true
#'   positive rates for each method.
#'
#' @examples
#'
#' pm = example1$pm
#' am = example1$am
#' missing = example1$missing
#' refs = "R1"
#'
#' db = forrel::NorwegianFrequencies[1:3]
#'
#' # True solution
#' true = c("M1", "M2", "M3")
#'
#' # Run comparison
#' \donttest{
#' dviCompare(pm, am, missing, refs, true = true, db = db, Nsim = 2, seed = 123)
#' }
#'
#' # Alternatively, simulations can be done first...
#' sims = dviCompare(pm, am, missing, refs, true = true, simulate = TRUE,
#'                   db = db, Nsim = 2, seed = 123, returnSims = TRUE)
#'
#' # ... and computations after:
#' \donttest{
#' dviCompare(sims$pm, sims$am, missing, refs, true = true, simulate = FALSE)
#' }
#' 
#' @importFrom forrel profileSim
#' @importFrom parallel makeCluster stopCluster parLapply clusterEvalQ
#'   clusterExport clusterSetRNGStream
#' @export
dviCompare = function(pm, am, missing, true, refs = typedMembers(am), methods = 1:6, 
                      markers = NULL, threshold = 1, simulate = TRUE, 
                      db = getFreqDatabase(am), Nsim = 1, returnSims = FALSE, 
                      seed = NULL, numCores = 1, verbose = TRUE) {
  st = Sys.time()
  
  if(is.singleton(pm))
    pm = list(pm)
  
  missing = as.character(missing)
  refs = as.character(refs)
  true = as.character(true)
  
  if(verbose) {
    if(simulate) 
      summariseDVI(pm, am, missing, printMax = 10)
    else
      summariseDVI(pm[[1]], am[[1]], missing, printMax = 10)
    message("\nParameters for DVI comparison:")
    message(" True solution: ", toString(true))
    message(" Simulate data: ", simulate)
    message(" Number of sims: ", if(simulate) Nsim else length(pm))
    message(" Reference IDs: ", toString(refs))
    message(" LR threshold: ", threshold)
    message("")
  }
  
  if(simulate) {
    vics = names(pm) = unlist(labels(pm))
    isMatch = true != "*"
    
    stopifnot(length(true) == length(vics), 
              all(true[isMatch] %in% missing), 
              all(getSex(am, true[isMatch]) == getSex(pm)[isMatch]))
    
    if(!is.null(seed))
      set.seed(seed)
    
    am = setMarkers(am, locusAttributes = db)
    pm = setMarkers(pm, locusAttributes = db)
    
    # Simulate AM
    AMsims = profileSim(am, N = Nsim, ids = c(refs, true[isMatch]))
    
    # Simulate the unrelated victims
    PMsims = forrel::profileSim(pm, N = Nsim, ids = vics[!isMatch])
    
    # For the true matches, transfer from MP to vics
    PMsims = lapply(1:Nsim, function(i) 
      transferMarkers(from = AMsims[[i]], to = PMsims[[i]], 
                      idsFrom = true[isMatch], idsTo = vics[isMatch], erase = FALSE))
    
    # Remove data from missing
    AMsims = lapply(AMsims, function(s) setAlleles(s, missing, alleles = 0))
    
    # Return sims
    if(returnSims) {
      return(list(pm = PMsims, am = AMsims))
    }
    
  }
  else {
    vics = unlist(labels(pm[[1]]))
    stopifnot(length(true) == length(vics), 
              all(true %in% c(missing, "*")),
              setequal(refs, typedMembers(am[[1]])))
    
    PMsims = pm
    AMsims = am
    if(!is.null(markers)) {
      PMsims = lapply(PMsims, selectMarkers, markers)
      AMsims = lapply(AMsims, selectMarkers, markers)
    }
  }
  
  N = length(PMsims)
  stopifnot(N == length(AMsims))
  
  # Setup parallelisation
  if(paral <- (numCores > 1)) {
    cl = makeCluster(numCores)
    if(verbose) 
      message("Using ", length(cl), " cores")
    on.exit(stopCluster(cl))
    clusterEvalQ(cl, library(dvir))
    clusterExport(cl, c("missing"), envir = environment())
    clusterSetRNGStream(cl, iseed = sample.int(1e6,1))  
  }
  
  # DVI functions (just to reduce typing)
  fun1 = function(i) sequentialDVI(PMsims[[i]], AMsims[[i]], missing, threshold = threshold, updateLR = FALSE, check = FALSE, verbose = FALSE)
  fun2 = function(i) sequentialDVI(PMsims[[i]], AMsims[[i]], missing, threshold = threshold, updateLR = TRUE, check = FALSE, verbose = FALSE)
  fun3 = function(i) top(jointDVI(PMsims[[i]], AMsims[[i]], missing, undisputed = TRUE, threshold = threshold, check = FALSE, verbose = FALSE))
  fun4 = function(i) top(jointDVI(PMsims[[i]], AMsims[[i]], missing, undisputed = FALSE, check = FALSE, verbose = FALSE))
  fun5 = function(i) top(jointDVI(PMsims[[i]], AMsims[[i]], missing, undisputed = TRUE, threshold = threshold, check = FALSE, verbose = FALSE), threshold)
  fun6 = function(i) top(jointDVI(PMsims[[i]], AMsims[[i]], missing, undisputed = FALSE, check = FALSE, verbose = FALSE), threshold)
  
  # Initialise list of results
  res = list()
  
  # Approach 1: Sequential - naive
  if(1 %in% methods) {
    method1 = if(paral) parLapply(cl, 1:N, fun1) else lapply(1:N, fun1)
    res$method1 = summar(method1)
    if(verbose) print(res['method1'])
  }
  
  # Approach 2: Sequential - with LR update
  if(2 %in% methods) {
    method2 = if(paral) parLapply(cl, 1:N, fun2) else lapply(1:N, fun2)
    res$method2 = summar(method2)
    if(verbose) print(res['method2'])
  }
  
  # Approach 3: Sequential undisputed + joint
  if(3 %in% methods) {
    method3 = if(paral) parLapply(cl, 1:N, fun3) else lapply(1:N, fun3); 
    res$method3 = summar(method3)
    if(verbose) print(res['method3'])
  }
  
  # Approach 4: Joint
  if(4 %in% methods) {
    method4 = if(paral) parLapply(cl, 1:N, fun4) else lapply(1:N, fun4)
    res$method4 = summar(method4)
    if(verbose) print(res['method4'])
  }
  
  # Approach 5: Sequential undisputed + joint + threshold (LR < thresh => *-*-*)
  if(5 %in% methods) {
    method5 = if(paral) parLapply(cl, 1:N, fun5) else lapply(1:N, fun5); 
    res$method5 = summar(method5)
    if(verbose) print(res['method5'])
  }
  
  # Approach 6: Joint  + threshold (LR < thresh => *-*-*)
  if(6 %in% methods) {
    method6 = if(paral) parLapply(cl, 1:N, fun6) else lapply(1:N, fun6)
    res$method6 = summar(method6)
    if(verbose) print(res['method6'])
  }
  
  # True positive rates
  TPR = sapply(res, function(r) as.numeric(r[paste(true, collapse = "-")]))
  TPR[is.na(TPR)] = 0
  res$TPR = TPR
  
  if(verbose)
    message("Total time used: ", format(Sys.time() - st, digits = 3))
  
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
  sort(freqs, decreasing = T)
}

# Utility for including ties in output of jointDVI()
top = function(res, thresh = 1) {
  nvic = ncol(res) - 3
  if(res$LR[1] < thresh) return(rbind(rep("*", nvic)))
  mx = res$LR == res$LR[1]
  res[mx, seq_len(nvic), drop = FALSE]
}
