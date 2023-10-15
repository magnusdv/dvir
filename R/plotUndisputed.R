#' Plot undisputed identifications
#' 
#' @inheritParams plotSolution
#' @param undisputed A data frame containing the undisputed matches, typically
#'   the entry `undisputed` in output from [findUndisputed()] (only three first
#'   columns used).
#' @param ... Further arguments passed on to [plotSolution()].
#'
#' @return NULL
#' 
#' @seealso [findUndisputed()], [plotSolution()]
#' @examples
#'
#' # Example
#' res = findUndisputed(example2, threshold = 2, verbose = FALSE)
#' u = res$summary
#' plotUndisputed(example2, u, marker = 1)
#'
#' @export
plotUndisputed = function(dvi, undisputed, ...){
  
  d = dim(undisputed)[1]
  if(d == 0)
    stop2("No undisputed to plot.")
  
  # Include only families with identification(s)
  dvi = subsetDVI(dvi, am = unique(undisputed$Family), verbose = FALSE)
  
  # Assignment as named vector (vic = miss)
  a = undisputed$Missing
  names(a) = undisputed$Sample

  plotSolution(dvi, assignment = a, ...)
}