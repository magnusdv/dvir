#' Simulate genotypes in a DVI dataset
#'
#' Simulates genotypes for the references and missing persons in each AM family,
#' transfers to the PM singletons according to the indicated matching. Remaining
#' victims are simulated as unrelated.
#'
#' @param dvi A `dviData` object.
#' @param refs A character with names of all reference individuals. By default,
#'   the typed members of the input.
#' @param truth A named vector of the format `c(vic1 = mis1, vic2 = mis2, ...)`.
#' @param seed An integer seed for the random number generator.
#' @param verbose A logical.
#'
#' @return A `dviData` object similar to the input.
#'
#' @seealso [forrel::profileSim()]
#'
#' @examples
#' ex = dviSim(example2, truth = c(V1 = "M1", V2 = "M2"))
#' plotDVI(ex, marker = 1)
#'
#' @export
dviSim = function(dvi, refs = typedMembers(dvi$am), truth = NULL, seed = NULL, 
                  verbose = FALSE){
  
  dvi = consolidateDVI(dvi)
  
  if(!is.null(seed))
    set.seed(seed)
  
  # Remove existing genotypes
  pm = dvi$pm |> setAlleles(alleles = 0)
  am = dvi$am |> setAlleles(alleles = 0)
  missing = dvi$missing
  
  # Simulate profiles of references and missing persons
  amSim = profileSim(am, ids = c(refs, missing))
  
  # Transfer to PM according to `truth`
  pmSim = transferMarkers(from = amSim, to = pm, idsFrom = truth,
                          idsTo = names(truth), erase = FALSE) 
  
  # Simulate the remaining victims
  remVics = setdiff(names(pm), names(truth))
  pmSim[remVics] = forrel::profileSim(pm[remVics])
  
  # Erase genotypes of missing
  amSim = setAlleles(amSim, ids = missing, alleles = 0)
  
  dviData(pmSim, amSim, missing)
}
