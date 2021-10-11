#' Age based DVI search
#'
#' Based on uncertain estimates of ages of victims, and certain age information 
#' on missing persons, the log likelihood, AIC and posterior 
#' are computed for all assignments-
#'
#' @param tab Data frame of assignments.
#' @param ageV List of victim ages.
#' @param ageM List of missing person ages.
#' @param k Double. Unconditional victim ages are assumed uniform on `[0,k]`.
#' @param sigma Conditional victim age V | V = M are assumed normal 
#'              with expectation equal to the missing person and standard 
#'              deviation sigma.
#' @param prior Double vector, default flat.             
#' @details The likelihood of victim ages conditioned on all assignments,
#' is calculated assuming assuming independence. The unconditional distribution 
#' of the age of victim ' (i.e., victim not identified) is uniform.
#' The conditional distribution is normal with expectation 
#' equal to the corresponding missing person and a user specified 
#' standard deviation sigma. 
#' The best model based on AIC is the model that minimises 
#' AIC = 2*|a| -2*log(likelihood) where |a| is the number of victims identified.
#'             
#'
#' @return A data frame. Each row describes an assignment of victims to missing
#'   persons, accompanied with its age based log likelihood, AIC and posterior
#' @seealso [expand.grid.nodup2()]
#'
#' @examples
#' pm = example2$pm
#' am = example2$am
#' missing = example2$missing
#' am[1] = setAlleles(am[1],"R1", alleles = 0)
#' pm[[3]] = swapSex(pm[[3]], "V3")
#' am[[2]] = swapSex(am[[2]], "M3")
#' pm = setMutationModel(pm, model = "proportional", rate = 0.01)
#' am = setMutationModel(am, model = "proportional", rate = 0.01)
#' miss = c('*', missing)
#' lst = list(V1 = miss, V2 = miss, V3 = miss)
#' tab = expand.grid.nodup2(lst, pm, am)
#' ageM = list(M1 = 58, M2 = 10, M3 = 20)
#' 
#' # Simulate ages from assignment (V1 = M1, V2 = M2)
#' set.seed(177)
#' sigma = 2
#' ageV = list(V1 = rnorm(1, 58, sigma), V2 = rnorm(1, 10, sigma), 
#'           V3 = runif(1, 0, 100))
#' res = ageAIC(tab, ageV, ageM, k = 100, sigma = 2)
#' # The assignment from which victim ages were simulated 
#' # comes out as the most likely
#' res[order(res$AIC, decreasing = F)[1:5],]
#' @export
#' 

ageAIC = function(tab, ageV, ageM, k = 100, sigma = 2, prior = NULL){
  # Add testing of input
  # Replace missing persons with their age
  age = convertToAge(tab, ageM)
  na = dim(tab)[1]
  nV = dim(tab)[2]
  # Flat prior if not otherwise specified
  if(is.null(prior))
    prior = rep(1/na, na)
  
  # Calculate for each line of assignment table:
  l = AIC =  rep(NA, na)
  for (i in 1:na){
    x = age[i,]
    star = x == '*'
    nstar = sum(star)
    aV = as.double(ageV[!star])
    aM = as.double(x[!star])
    l[i] = -nstar*log(k)+sum(dnorm(aV, aM, sigma, log = T))
    AIC[i] = 2*(nV-nstar)-2*l[i]
  }
  lmax = max(l)
  term = prior * exp(l- lmax)
  posterior = term/sum(term)
  data.frame(tab, prior, logLikAge = l,  posterior, AIC)
}



# Function to replace missing person by their age
convertToAge = function(tab, ageM){
  age = tab
  na = dim(tab)[1] 
  nV = dim(tab)[2] 
  for (i in 1:na)
    for (j in 1:nV)
      if (tab[i,j] != "*")
        age[i,j] = ageM[tab[i,j]]
  age
}

# @rdname ageAIC
# @export

