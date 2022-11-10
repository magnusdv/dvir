#' @importFrom pedprobr likelihood
loglikTotal = function(x, markers = seq_len(nMarkers(x))) {
  sum(likelihood(x, marker = markers, logbase = exp(1), eliminate = 1))
}


#' @importFrom pedprobr setMutationModel
#' @export
pedprobr::setMutationModel

# Modified version of stop()
stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}


`%||%` = function(x, y) {
  if(is.null(x)) y else x
}
