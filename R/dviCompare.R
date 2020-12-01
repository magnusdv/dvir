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
#' @param db A frequency database, e.g., forrel::NorwegianFrequencies[1:10].
#' @param Nsim A positive integer.
#' @param methods A subset of the numbers 1,2,3.
#' @param seed A seed for the random number generator, or NULL.
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
#' dviCompare(pm, am, MPs, refs, 
#'            true = true, db = db, Nsim = 2, seed = 1234)
#' 
#' @importFrom forrel profileSim
#' @export
dviCompare = function(pm, am, MPs, refs, true, db, Nsim, methods = 1:3, seed = NULL, verbose = FALSE) {

  vics = names(pm) = unlist(labels(pm))
  
  am = setMarkers(am, locusAttributes = db)
  pm = setMarkers(pm, locusAttributes = db)
  
  # Which victims are actually MPs
  isMatch = true != "*"
  
  if(!is.null(seed))
    set.seed(seed)
  
  # Simulate AM
  AMsims = profileSim(am, N = Nsim, ids = c("R1", MPs))
  
  # Simulate the unrelated victims
  PMsims = profileSim(pm, N = Nsim, ids = vics[!isMatch])
  
  # For the true matches, transfer from MP to vics
  PMsims = lapply(1:Nsim, function(i) 
    transferMarkers(from = AMsims[[i]], to = PMsims[[i]], 
                    idsFrom = true[isMatch], idsTo = vics[isMatch]))
  
  # Remove data from missing
  AMsims = lapply(AMsims, function(s) setAlleles(s, MPs, alleles = 0))
  
  # Initialise list of results
  res = list()
  
  # Approach 1: Sequential - naive
  if(1 %in% methods) {
    seq1 = lapply(1:Nsim, function(i) sequential1(PMsims[[i]], AMsims[[i]], MPs, verbose = FALSE))
    res$seq1 = summar(seq1)
    if(verbose) print(res['seq1'])
  }
  
  # Approach 2: Sequential - with replacement
  if(2 %in% methods) {
    seq2 = lapply(1:Nsim, function(i) sequential2(PMsims[[i]], AMsims[[i]], MPs, verbose = FALSE))
    res$seq2 = summar(seq2)
    if(verbose) print(res['seq2'])
  }
    # Approach 3: Joint
  if(3 %in% methods) {
    joint = lapply(1:Nsim, function(i) 
      global(PMsims[[i]], AMsims[[i]], MPs, check = FALSE, verbose = F)[1, 1:length(vics)])
    res$joint = summar(joint)
    if(verbose) print(res['joint'])
  }
  
  # True positive rates
  TPR = sapply(res, function(r) as.numeric(r[paste(true, collapse = "-")]))
  TPR[is.na(TPR)] = 0
  
  res$TPR = TPR
  
  res
}

# Utility for summarising list of solutions
summar = function(x) {
  tab = table(sapply(x, function(y) paste(y, collapse = "-")))
  freqs = as.vector(tab/length(x))
  names(freqs) = names(tab)
  sort(freqs, decreasing = T)
}

