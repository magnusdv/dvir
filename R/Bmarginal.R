#' Posterior pairing probabilities
#'
#' Compute posterior pairing and non-pairing probabilities, based on a prior and
#' the output from [jointDVI()].
#'
#' The prior assigns a probability to each assignment, each row of `jointRes`.
#' If the prior is not specified, a flat prior is used. The prior needs not sum
#' to 1 since the user may rather choose a flat prior on the *a priori* possible
#' assignments.
#'
#' @param jointRes Output from [jointDVI()].
#' @param missing Character vector with names of missing persons.
#' @param prior A numeric vector of length equal the number of rows in
#'   `jointRes`. Default is a flat prior.
#'
#' @return A matrix. Row `i` gives the posterior probability that victim `i` is
#'   one of the missing persons or someone else, denoted '*'.
#'
#' @seealso [jointDVI()]
#'
#' @examples
#' jointRes = jointDVI(example1)
#'
#' Bmarginal(jointRes, example1$missing)
#'
#' # Artificial example: all but optimal solution excluded by prior
#' Bmarginal(jointRes, example1$missing, prior = c(1, rep(0,26)))
#'
#'
#' @export
Bmarginal = function(jointRes, missing, prior = NULL){
  
  d = dim(jointRes)
  
  # Deal with prior
  if(is.null(prior))
    prior = rep(1/d[1], d[1])  
  else if(length(prior) != d[1])
    stop2("Length of prior must equal number of rows in first argument")
  
  # Find names of victims
  nv = (1:d[2])[names(jointRes) == "loglik"] - 1
  victims = names(jointRes)[1:nv]
  
  # Initialise result table 
  val = c(missing, '*')
  res = matrix(nrow = nv, ncol = length(val), 0)
  dimnames(res) = list(victims, val)
  
  # Find denominator, const, and do calculations
  # Subtract max loglik in numerator and denominator to avoid underflow
  lmax = max(jointRes$loglik) 
  x = data.frame(jointRes, term = prior * exp(jointRes$loglik-lmax))
  const = sum(x$term)  
  for(v in 1:nv){
    r1 = split(x$term, x[,v])
    p = sapply(r1, sum)/const
    res[v, names(p)] = p
  }
  
  res
}