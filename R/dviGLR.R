#' Finds Generalised Likelihood Ratios (GLRs)
#'
#' Based on a `dviData` object, or output from `dviJoint()``or `jointDVI()`, 
#' the GLR, the ratio of the maximum likelihood under H0 to the maximum under H1,
#' is calculated for specified hypotheses.
#'
#' @param dvi A `dviData` object.
#' @param pairings List. See details.
#' @param dviRes data frame. Output from `jointDVI()` or `dviJoint()`.
#'
#' @return A data frame with GLRs and SGLR (strict GLR,max replaced by min in 
#'   the numerator).
#'
#' @details The Generalised Likelihood Ratio (GLR) statistic is defined as the
#'   ratio of the maximum likelihood for the alternatives in the numerator
#'   to the maximum in the denominator. The default 
#'   `pairings = dvir::generatePairings(dvi)` tests all hypotheses. Specific 
#'   tests can be specified as shown in an example:
#'   `pairings = list(V1 = "M1")` gives a test for H0: V1 = M1 against 
#'   H1: V1 != M1. `dviRes` will be calculated using `dviJoint()` if not provided.
#'
#' @examples
#' dviGLR(example2, pairings = list(V1 = "M1"))
#'
#' r = jointDVI(example2)
#' dviGLR(example2, pairings = list(V1 = "M1"), dviRes = r)
#'
#' # All tests with output from dviJoint
#' dviGLR(example2, dviRes = r)
#'
#' @export
dviGLR = function(dvi, pairings = generatePairings(dvi), dviRes = NULL){
  if(!inherits(dvi, "dviData"))
    stop("First argument needs to be a dviData object")

  # Do joint analysis, if not provided, and find max loglikelihood
  if(is.null(dviRes))
    dviRes = jointDVI(dvi)
  loglik = dviRes$loglik
  na = dim(dviRes)[1]

  # The number of hypotheses to be tested
  nhyp = sum(unlist(lapply(pairings, function(x) length(x))))

  # Initialise log likelihoods, max and min under H0, and max under H1
  l0Max = l0Min = l1 = rep(-Inf, nhyp)

  # To contain names of null hypotheses in output:
  hyp = rep(NA, nhyp)

  # Loop through hypotheses for each victim
  npairs = length(pairings)
  s = 0
  for (i in 1:npairs){
    victim = names(pairings[i])
    for(j in pairings[[i]]){
      s = s + 1
      hyp[s] = paste(victim, j, sep = "=")
      
      # Find indices for numerator loglikelihoods
      I0 = (1:na)[dviRes[,victim] == j]
      I1 = setdiff(1:na, I0)
      if(length(I0)) {
        l0Max[s] = max(loglik[I0])
        l0Min[s] = min(loglik[I0])
      }
      if(length(I1)) l1[s] = max(loglik[I1])
    }
  }
  GLR = exp(l0Max - l1)
  SGLR = exp(l0Min - l1)
  data.frame(H0 = hyp, GLR = GLR, GLRmin = SGLR)
}
