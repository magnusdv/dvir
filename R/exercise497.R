#' Solution Exercise 4.9.7 in the book Kling et al. (2021)
#'
#' This is a DVI case with 3 female victims and 3 missing females
#' in three reference families. There are 23 markers with equal mutation rate 0.001.
#' Data are simulated from the solution V1 = MP1, V2 = MP2, V3 = MP3 and the purpose is
#' to check fraction of times the 'correct' solutions is obtained.
#'
#' @param pm A list of singletons.
#' @param am A list of pedigrees.
#' @param missing Character vector with names of missing persons.
#' @param nsim Number of simulations.
#' @param seed Integer.
#' @param simRef Logical. if TRUE, references are also simulated.
#' @param disableMutations Logical, see [jointDVI()].
#' @param undisputed Logical, see [jointDVI()].
#' @param verbose A logical.
#
#'
#' @return A list with two elements, the first the fraction of 'correct' solutions, the second
#' a matrix with first line from [jointDVI()]
#'
#' @seealso [jointDVI()]
#'
#' @examples
#' \donttest{
#' pm = dataExercise497$pm
#' am = dataExercise497$am
#' missing = dataExercise497$missing
#' exercise497(pm, am, missing, nsim = 2, seed = 17, verbose = TRUE)
#' }
#' @export
#' 
exercise497 = function(pm, am, missing, nsim = 2, seed = NULL, 
                       simRef = TRUE,  disableMutations = FALSE, 
                       undisputed = FALSE, verbose = FALSE){
  set.seed(seed)
  
  # Simulate V1 = MP1, V2 = MP2, V3 = MP3, and also references if simRef = TRUE 
  if (simRef){
    typed = typedMembers(am)
    am2 = setAlleles(am, alleles = 0)
    amsim = profileSim(am2, N = nsim , c(missing, typed))
  } else 
    amsim = profileSim(am, N = nsim , missing)
  
  victims = paste0("V", 1:3)
  resAll = NULL
  ncorrect = 0
  
  for (i in 1:nsim){
    pmsim = transferMarkers(amsim[[i]], pm, idsFrom = missing, idsTo = victims)
    am0 = setAlleles(amsim[[i]], ids = missing, alleles = 0)
    res = jointDVI(pmsim, am0, missing, disableMutations = disableMutations, 
                   undisputed = undisputed, verbose = FALSE)
    resAll = rbind(resAll, res[1,])
    if (verbose) {
      cat("Simulation: ", i, "\n")
      print(res[1,])
    }
    ncorrect = ncorrect + all(as.character(res[1,1:3]) == missing)
  }
  list(fractionCorrect = ncorrect/nsim, resAll =  resAll)
}
