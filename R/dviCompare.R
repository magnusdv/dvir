#' Compare DVI approaches
#'
#' The following methods are compared: 1 = sequential (simple); 2 = sequential
#' (with updates); 3 = joint.
#'
#' @param pm PM data: List of singletons
#' @param am AM data: A ped object or list of such.
#' @param MPs Character vector with names of the missing persons.
#' @param refs Character vector with names of the reference individuals.
#' @param true A character of the same length as `pm`, with the true solution,
#'   e.g., `true = c("MP2", "*", "MP3)` if the truth is V1 = MP2 and V3 = MP3.
#' @param methods A subset of the numbers 1,2,3.
#' @param markers If `simulate = FALSE`: A vector indicating which markers
#'   should be used.
#' @param simulate A logical, indicating if simulations should be performed.
#' @param db A frequency database used for simulation, e.g.,
#'   forrel::NorwegianFrequencies[1:10].
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
#' @importFrom parallel makeCluster stopCluster parLapply clusterEvalQ clusterExport clusterSetRNGStream
#' @export
dviCompare = function(pm, am, MPs, refs, true, methods = 1:3, markers = NULL, 
                      simulate = TRUE, db = NULL, Nsim = 1, returnSims = FALSE, 
                      seed = NULL, numCores = 1, verbose = FALSE) {
  st = Sys.time()
  
  if(simulate) {
    vics = names(pm) = unlist(labels(pm))
    isMatch = true != "*"
    
    stopifnot(length(true) == length(vics), all(true[isMatch] %in% MPs))
    
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
    stopifnot(length(true) == length(vics), all(true %in% c(MPs, "*")),
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
  seq1Fun = function(i) sequential1(PMsims[[i]], AMsims[[i]], MPs, check = FALSE, verbose = FALSE)
  seq2Fun = function(i) sequential2(PMsims[[i]], AMsims[[i]], MPs, check = FALSE, verbose = FALSE)
  jointFun = function(i) global(PMsims[[i]], AMsims[[i]], MPs, check = FALSE, verbose = F)[1, 1:length(vics)]
  
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
    # Approach 3: Joint
  if(3 %in% methods) {
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

