#' Simulate and evaluate DVI cases
#'
#' Individuals in the am data are first simulated. 
#' The marker data for the missing persons are then transferred to the pm data according
#' to the user specified solution to the identification problem. Remaining victims are simulated
#' as unrelated. The fraction of times the 'correct'
#' (the one from which data was simulated) is estimated.
#'
#' @param pm A list of singletons.
#' @param am A list of pedigrees.
#' @param missing Character vector with names of missing persons.
#' @param simulated Character vector with names of persons in am to be simulated
#' @param true A character of the same length as `pm`, with the true solution,
#' e.g., `true = c("M2", "*", "M3)` if the truth is V1 = M2 and V3 = M3.
#' @param Nsim Number of simulations.
#' @param seed Integer.
#' @param disableMutations Logical, see [jointDVI()].
#' @param undisputed Logical, see [jointDVI()].
#' @param verbose A logical
#
#' @return A list with two elements, the first the fraction of 'correct' solutions, the second
#' a matrix with first line, best solution, from [jointDVI()]
#'
#' @seealso [jointDVI()]
#'
#' @examples
#' \donttest{
#' library(pedtools)
#' pm = dataExercise497$pm
#' pm = setAlleles(pm, alleles = 0)
#' am = dataExercise497$am
#' missing = dataExercise497$missing
#' true = c("MP1", "MP2", "MP3")
#' true = c("MP1", "*", "MP3")
#' true = c("*", "*", "*")
#' dviSim(pm, am, missing, true = true, Nsim = 2, seed = 17, verbose = TRUE)
#' 
#' pm = example2$pm
#' pm = setAlleles(pm, alleles = 0)
#' am = example2$am
#' missing = example2$missing
#' true = c("M1", "*", "*")
#' true = c("M1", "M2", "*")
#' true = c("M1", "M2", "M3")
#' dviSim(pm, am, missing, true = true, Nsim = 2, seed = 17, verbose = FALSE)
#' 
#' pm = example1$pm
#' pm = setAlleles(pm, alleles = 0)
#' am = example1$am
#' missing = example1$missing
#' true = c("M1", "*", "*")
#' true = c("M1", "M2", "*")
#' true = c("M1", "M2", "M3")
#' dviSim(pm, am, missing, true = true, Nsim = 10, seed = 17, verbose = FALSE)
#' }
#' @export
dviSim = function(pm, am, missing, simulated = missing, true = NULL, Nsim = 2, seed = NULL, 
                  disableMutations = NA, undisputed = TRUE, verbose = FALSE){
  
  # Get names of victims and check consistency of 'true' assignment
  vics = names(pm) = unlist(labels(pm))
  isMatch = true != "*"
  stopifnot(length(true) == length(vics), 
            all(true[isMatch] %in% missing), 
            all(getSex(am, true[isMatch]) == getSex(pm)[isMatch]))
  
  # Initialise variables to contain results
  resAll = NULL
  ncorrect = 0
  
  # Simulate MP-s and possibly references
  amsim = profileSim(am, N = Nsim, simulated, seed = seed)
  
  # Simulate the unrelated victims
  pmsim = forrel::profileSim(pm, N = Nsim, ids = vics[!isMatch])
  
  for (i in 1:Nsim){
    # For the true matches, transfer from MP to vics
    pmsim[[i]] = transferMarkers(amsim[[i]], pmsim[[i]], 
                            idsFrom = true[isMatch], idsTo = vics[isMatch], erase = FALSE)
    am0 = setAlleles(amsim[[i]], ids = missing, alleles = 0)
    res = jointDVI(pmsim[[i]], am0, missing, disableMutations = disableMutations, 
                   undisputed = undisputed, verbose = FALSE)
    resAll = rbind(resAll, res[1,])
    if (verbose) {
      cat("Simulation: ", i, "\n")
      print(res[1,])
    }
    ncorrect = ncorrect + all(as.character(res[1,1:3]) == true)
  }
  list(fractionCorrect = ncorrect/Nsim, resAll =  resAll)
}
