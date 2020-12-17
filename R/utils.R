#' @importFrom pedprobr likelihood
loglikTotal = function(x, markers = seq_len(nMarkers(x))) {
  sum(likelihood(x, marker = markers, logbase = exp(1), eliminate = 1))
}


#' @importFrom pedprobr setMutationModel
#' @export
pedprobr::setMutationModel
