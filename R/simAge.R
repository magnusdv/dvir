#' Simulation experiment to evaluate age information
#'
#' Missing person ages are simulated, default uniformly on `(0, k)`,
#' but obeying pedigree ordering.
#' Victim ages are simulated conditional on assignment. If a victim 
#' is not identified by assignment, age is uniform on `(0, k)`.
#' Otherwise, age of V given V = M, is normal with expectation equal to
#' age of M and standard deviation sigma. The purpose is to see how well
#' a given assignment is found based on likelihood and AIC. Also, the posterior is
#' reported and can be compared to the prior.
#'
#' @param tab Data frame of assignments.
#' @param missing Character vector giveing names of the missing persons.
#' @param ageM List with ages of missing persons.
#' @param k Double. Unconditional victim ages are assumed uniform on `[0,k]`.
#' @param sigma Double.
#' @param prior Double vector, default `NULL` gives flat prior.
#'             
#'
#' @return A data frame. Each row describes an assignment of victims to missing
#'   persons, accompanied with its age based log likelihood, AIC, posterior and a column correct
#'   indicating with a 1 if the assignment from which simulation was done is found.
#' @seealso [ageAIC()]
#'
#' @examples
#' # Load example
#' pm = example2$pm
#' am = example2$am
#' missing = example2$missing
#' am[1] = setAlleles(am[1],"R1", alleles = 0)
#' pm[[3]] = swapSex(pm[[3]], "V3")
#' am[[2]] = swapSex(am[[2]], "M3")
#' # Generate assignments
#' miss = c('*', missing)
#' lst = list(V1 = miss, V2 = miss, V3 = miss)
#' tab = expand.grid.nodup2(lst, pm, am)
#' M1 = runif(1, 20, 100)
#' sigma = 2
#' M2 = max(M1-20-rnorm(1,0,sigma), 1)
#' M3 = runif(1, 20, 80)
#' ageM = list(M1 = M1, M2 = M2, M3 = M3)
#' res = simAge(tab, missing, ageM)
#' # Fraction correctly identified:
#' correct = res$correct
#' sum(correct)/length(correct)
#' \dontrun{
#' date()
#' nsim = 100
#' fracCorrect = rep(nsim, nsim)
#' set.seed(1729)
#' res = list()
#' for (i in 1:nsim){
#'   M1 = runif(1, 20, 100)
#'   M2 = max(M1-20-rnorm(1, 0, sigma), 1)
#'   M3 = runif(1, 20, 80)
#'   ageM = list(M1 = M1, M2 = M2, M3 = M3)
#'   res[[i]] = simAge(tab, missing, ageM)
#'   correct = res[[i]]$correct
#'   fracCorrect[i] = sum(correct)/length(correct)
#' }
#' fracCorrect
#' date()
#' }

#' @export
#' 

simAge = function(tab, missing, ageM, k = 100, sigma = 2, prior = NULL){
  na = dim(tab)[1]
  posterior = logLik = AIC = correct = rep(NA, na)
  if(is.null(prior))
    prior = rep(1/na, na)
  for (i in 1:na){
    a = tab[i,]
    aV = simV(a, ageM, sigma = sigma)
    res = ageAIC(tab, aV, ageM, k = k, sigma = sigma)
    logLik[i] = res$logLikAge[i]
    AIC[i] = res$AIC[i]
    correct[i] = (which.min(res$AIC) == i)
    posterior[i] = res$posterior[i]
  }
  data.frame(tab, prior, logLik, posterior, AIC,correct)
}


# Simulate victim ages for an assignment
simV = function(a, ageM, k = 100, sigma = 2){
  nV = length(a)
  b = as.list(a)
  aV = rep(NA, nV)
  for (i in 1:nV){
    if(a[i] == '*')
      aV[i] = runif(1, 0, k)
    else
      aV[i] = rnorm(1,ageM[[b[[i]]]], sigma)
  }
  aV
}


