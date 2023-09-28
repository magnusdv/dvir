#' Plot undisputed identifications
#'
#' @param dvi A `dviData` object, typically created with [dviData()].
#' @param undisputed A data frame containing the undisputed matches, typically
#'   the entry `undisputed` in output from [findUndisputed()] 
#'   (only three first columns used).
#' @param format See [plotSolution()]. 
#' @param ... Further parameters to be passed on to [pedtools::plot.ped()],
#'   e.g., `marker`, `cex`, `cex.main`, `symbolsize`.
#'
#' @return NULL
#'
#' @examples
#'
#' # Example 
#' res = findUndisputed(example2, threshold = 2, relax = TRUE, verbose = FALSE)
#' u = res$undisputed
#' plotUndisputed(example2, u, marker = 1)
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics grconvertX grconvertY layout mtext par rect
#' @export
plotUndisputed = function(dvi, undisputed,  format = "[S]=[M]", ...){
  
  d = dim(undisputed)[1]
  if(d == 0)
    stop("No undisputed to plot.")
  
  #Include only families with identification(s)
  dvi$am = dvi$am[unique(undisputed[,3])]
  
  x = as.data.frame(matrix(nrow = 1, ncol = d))
  x[1,] = undisputed[,2]
  colnames(x) = undisputed[,1]
  plotSolution(dvi, x, format = format, ...)
}