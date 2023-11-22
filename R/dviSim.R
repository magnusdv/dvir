#' Simulate genotypes in a DVI dataset
#'
#' Simulates genotypes for the references and missing persons in each AM family,
#' transfers to the PM singletons according to the indicated matching. Remaining
#' victims are simulated as unrelated.
#'
#' @param dvi A `dviData` object.
#' @param N	The number of complete simulations to be performed.
#' @param refs A character with names of all reference individuals. By default,
#'   the typed members of the input.
#' @param truth A named vector of the format `c(vic1 = mis1, vic2 = mis2, ...)`.
#' @param seed An integer seed for the random number generator.
#' @param conditional A logical, by default FALSE. If TRUE, references are kept
#'   unchanged, while the missing persons are simulated conditional on these.
#' @param simplify1	A logical, by default TRUE, removing the outer list layer 
#'   when N = 1. See Value.
#' @param verbose A logical.
#'
#' @return If `N = 1`, a `dviData` object similar to the input, but with new genotypes
#'   for the pm samples. If `N > 1`, a list of `dviData` objects. 
#'
#' @seealso [forrel::profileSim()].
#'
#' @examples
#'
#' # Simulate refs and missing once and plot:
#' ex = dviSim(example2, N = 1, truth = c(V1 = "M1", V2 = "M2", V3 = "M3"))
#' plotDVI(ex, marker = 1)
#'
#' # Two simulations and plot for the first
#' ex = dviSim(example2, N = 2, truth = c(V1 = "M1", V2 = "M2", V3 = "M3"), 
#'             seed = 1729)
#' plotDVI(ex[[1]], marker = 1)
#' 
#'
#' @export
dviSim = function(dvi, N = 1, refs = typedMembers(dvi$am), truth = NULL, 
                  seed = NULL, conditional = FALSE, simplify1 = TRUE,
                  verbose = FALSE) {
  dvi = consolidateDVI(dvi)
  
  labsAM = labels(dvi$am) |> unlist(use.names = FALSE)
  if (!all(refs %in% labsAM)) {
    stop2("Unknown reference ID of reference: ", setdiff(refs, labsAM))
  }
  
  if (any(refs %in% dvi$missing)) {
    stop2("Missing person cannot be a reference: ", intersect(refs, dvi$missing))
  }
  
  if (!all(truth %in% dvi$missing)) {
    stop2("Unknown missing person ID: ", setdiff(truth, dvi$missing))
  }
  
  if (!all(names(truth) %in% names(dvi$pm))) {
    stop2("Unknown victim ID: ", setdiff(names(truth), names(dvi$pm)))
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Remove existing genotypes in pm
  pm = dvi$pm |> setAlleles(alleles = 0)
  missing = dvi$missing
  
  # Simulate profiles of missing persons and also of references if conditional = FALSE
  if (conditional) {
    am = dvi$am |> setAlleles(alleles = 0, ids = missing)
    amSim = profileSim(am, N = N, ids = missing, simplify1 = FALSE, 
                       verbose = verbose)
  } else {
    am = dvi$am |> setAlleles(alleles = 0)
    amSim = profileSim(am, N = N, ids = c(refs, missing), simplify1 = FALSE, 
                       verbose = verbose)
  }
  
  # Transfer to PM according to `truth`
  pmSim = lapply(amSim, function(z) {
    transferMarkers(
      from = z, to = pm, idsFrom = truth, idsTo = names(truth), erase = FALSE
      )})
  
  # Erase genotypes of missing
  amSim = lapply(amSim, function(z) setAlleles(z, ids = missing, alleles = 0))
  
  # Simulate the remaining victims
  remVics = setdiff(names(pm), names(truth))
  if (length(remVics)) 
    pmSim = lapply(pmSim, function(z) forrel::profileSim(z, ids = remVics))
  
  # Collect the list of dvi data objects  
  dviDataSimulated = lapply(seq_len(N), function(i) {
    dviData(pmSim[[i]], amSim[[i]], missing)
    })
 
  if (simplify1 && N == 1) 
    dviDataSimulated = dviDataSimulated[[1]]
  
  dviDataSimulated
}
