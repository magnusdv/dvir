#' Compare DVI approaches
#'
#' The following methods are compared: 1 = sequential (simple); 2 = sequential
#' (with updates); 3 = sequential + joint; 4 = joint.
#'
#' @param pm PM data: List of singletons
#' @param am AM data: A ped object or list of such.
#' @param MPs Character vector with names of the missing persons.
#' @param true A character of the same length as `pm`, with the true solution,
#'   e.g., `true = c("MP2", "*", "MP3)` if the truth is V1 = MP2 and V3 = MP3.
#' @param refs Character vector with names of the reference individuals. By
#'   default the typed members of `am`.
#' @param methods A subset of the numbers 1,2,3,4.
#' @param markers If `simulate = FALSE`: A vector indicating which markers
#'   should be used.
#' @param threshold An LR threshold passed on to the sequential methods.
#' @param simulate A logical, indicating if simulations should be performed.
#' @param db A frequency database used for simulation, e.g.,
#'   forrel::NorwegianFrequencies. By default the frequencies attached to `am`
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
#' MPs = example1$MPs
#' refs = "R1"
#'
#' db = forrel::NorwegianFrequencies[1:10]
#'
#' # True solution
#' true = c("MP1", "MP2", "MP3")
#'
#' # Run comparison
#' dviCompare(pm, am, MPs, refs, true = true,
#'            db = db, Nsim = 2, seed = 123)
#'
#'
#' # Alternatively, simulations can be done first...
#' sims = dviCompare(pm, am, MPs, refs, true = true,
#'                   db = db, Nsim = 2, seed = 123, returnSims = TRUE)
#'
#'  # ... and computations after:
#' dviCompare(sims$pm, sims$am, MPs, refs, true = true, simulate = FALSE)
#'
#' @importFrom forrel profileSim
#' @importFrom parallel makeCluster stopCluster parLapply clusterEvalQ
#'   clusterExport clusterSetRNGStream
#' @export
dviCompare = function(pm, am, MPs, true, refs = typedMembers(am), methods = 1:4, 
                      markers = NULL, threshold = 1, simulate = TRUE, 
                      db = getFreqDatabase(am), Nsim = 1, returnSims = FALSE, 
                      seed = NULL, numCores = 1, verbose = FALSE) {
  st = Sys.time()
  
  if(is.singleton(pm))
    pm = list(pm)
  
  MPs = as.character(MPs)
  refs = as.character(refs)
  true = as.character(true)
                      
  if(simulate) {
    vics = names(pm) = unlist(labels(pm))
    isMatch = true != "*"
    
    stopifnot(length(true) == length(vics), 
              all(true[isMatch] %in% MPs), 
              all(getSex(am, true[isMatch]) == getSex(pm)[isMatch]))
    
    if(!is.null(seed))
      set.seed(seed)
    
    am = setMarkers(am, locusAttributes = db)
    pm = setMarkers(pm, locusAttributes = db)
    
    # Simulate AM
    AMsims = profileSim(am, N = Nsim, ids = c(refs, true[isMatch]))
    
    # Simulate the unrelated victims
    PMsims = profileSim(pm, N = Nsim, ids = vics[!isMatch])
    
    # For the true matches, transfer from MP to vics
    PMsims = lapply(1:Nsim, function(i) 
      transferMarkers(from = AMsims[[i]], to = PMsims[[i]], 
                      idsFrom = true[isMatch], idsTo = vics[isMatch], erase = FALSE))
    
    # Remove data from missing
    AMsims = lapply(AMsims, function(s) setAlleles(s, MPs, alleles = 0))
    
    # Return sims
    if(returnSims) {
      return(list(pm = PMsims, am = AMsims))
    }
    
  }
  else {
    vics = unlist(labels(pm[[1]]))
    stopifnot(length(true) == length(vics), 
              all(true %in% c(MPs, "*")),
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
    clusterExport(cl, c("MPs"), envir = environment())
    clusterSetRNGStream(cl, iseed = sample.int(1e6,1))  
  }
  
  # DVI functions (just to reduce typing)
  seq1Fun = function(i) sequential1(PMsims[[i]], AMsims[[i]], MPs, threshold = threshold, check = FALSE)
  seq2Fun = function(i) sequential2(PMsims[[i]], AMsims[[i]], MPs, threshold = threshold, check = FALSE)
  seq3Fun = function(i) pickWinner(sequential3(PMsims[[i]], AMsims[[i]], MPs, threshold = threshold, check = FALSE))
  jointFun = function(i) pickWinner(global(PMsims[[i]], AMsims[[i]], MPs, check = FALSE))
  
  # Initialise list of results
  res = list()
  
  # Approach 1: Sequential - naive
  if(1 %in% methods) {
    seq1 = if(paral) parLapply(cl, 1:N, seq1Fun) else lapply(1:N, seq1Fun)
    res$seq1 = summar(seq1)
    if(verbose) print(res['seq1'])
  }
  
  # Approach 2: Sequential - with replacement
  if(2 %in% methods) {
    seq2 = if(paral) parLapply(cl, 1:N, seq2Fun) else lapply(1:N, seq2Fun)
    res$seq2 = summar(seq2)
    if(verbose) print(res['seq2'])
  }
  
  # Approach 3: Sequential undisputed + joint
  if(3 %in% methods) {
    seq3 = if(paral) parLapply(cl, 1:N, seq3Fun) else lapply(1:N, seq3Fun)
    res$seq3 = summar(seq3)
    if(verbose) print(res['seq3'])
  }
  
  # Approach 3: Joint
  if(4 %in% methods) {
    joint = if(paral) parLapply(cl, 1:N, jointFun) else lapply(1:N, jointFun)
    res$joint = summar(joint)
    if(verbose) print(res['joint'])
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
summar = function(x) {
  tab = table(sapply(x, function(y) paste(y, collapse = "-")))
  freqs = as.vector(tab/length(x))
  names(freqs) = names(tab)
  sort(freqs, decreasing = T)
}

# Utility for breaking ties in output of global()
pickWinner = function(res) {
  mx = which(res$LR == res$LR[1])
  
  # If winner is not unique, pick random 
  if(length(mx) > 1)
    mx = sample(mx, size = 1)
  
  # Return best assignment
  res[mx, seq_len(ncol(res) - 3)]
}