#' Marginal Bayesian identification probabilities
#'
#' Based on a prior and  the output from [jointDVI()], the posterior identification 
#' probabilities are found.
#'
#' @param jointRes Output from [jointDVI()].
#' @param missing Character vector with names of missing persons.
#' @param prior A numeric vector of length equal the number of rows on `jointRes`. 
#' Default is a flat prior.
#'
#' @return A matrix. Row i gives the posterior probability that victim i is 
#' one of the missimg persons or someone else, denoted '*'.
#' 
#' @details The prior assigns a probability to each assignment, each row of `jointRes`. 
#' If the prior is not specified, a flat prior is used. The prior need not sum to 1 since the user may
#' rather choose a flat prior on the `apriori` possible assignments. 
#'
#' @seealso [jointDVI()]
#'
#' @examples
#' data(example1)
#' pm = example1$pm
#' am = example1$am
#' missing = example1$missing
#' jointRes = jointDVI(pm, am, missing)
#' Bmarginal(jointRes,  missing)
#' Artificial example, all but optimal solution excluded by prior:
#' Bmarginal(jointRes,  missing, prior = c(1, rep(0,26)))
#' 
#' # Another example:
#' data(planecrash)
#' pm = planecrash$pm
#' am = planecrash$am
#' missing = planecrash$missing
#' jointRes = jointDVI(pm, am, missing)
#' Bmarginal(jointRes,  missing)
#' 
#'
#' @importFrom parallel makeCluster stopCluster detectCores parLapply
#'   clusterEvalQ clusterExport
#'
#' @export
Bmarginal = function(jointRes,  missing, prior = NULL){
  # Deal with prior
  if(!is.null(prior) & length(prior) != d[1])
    return("Length of prior should equal number of rows in first argument")
  if(is.null(prior))
    prior = rep(1/d[1], d[1])  
  # Find names of victims
  d = dim(jointRes)
  nv = (1:d[2])[names(jointRes)=="loglik"]-1
  victims = names(jointRes)[1:nv]
  # Initialise result table 
  val = c(missing, '*')
  res = matrix(nrow = nv, ncol = length(val), 0)
  dimnames(res) = list(victims, val)
  # Find denominator, const, and do calculations
  x = data.frame(jointRes, term = prior*exp(jointRes$loglik))
  const = sum(x$term)  
  for(v in 1:nv){
    r1 = split(x$term, x[,v])
    p = sapply(r1,sum)/const
    res[v, names(p)] = p
  }
  res
}