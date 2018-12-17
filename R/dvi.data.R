#' Data for dvir2: Simulated
#' 
#'  A list with 4 elements: "pm" (post mortem), "vict" (names of victims),
#' "am" (ante mortem), "miss" (names of missing persons). "pm" and "am" are ped objects
#' or lists of such while "vict" and "miss" are character vectors giving names of victims
#' and missing persons.

#' @name dvi.data
#' @docType data
#' @keywords data
#' @examples 
#' 
#' data(dvi.data)
#' pm = dvi.data$pm
#' am = dvi.data$am
#' library(pedtools)
#' plotPedList(list(pm,am), newdev = TRUE, marker = 1, 
#' skip.empty.genotypes = TRUE, dev.width = 10, dev.height = 7, cex =0.7, 
#' frametitles = c("Post mortem data (first marker)", "Ante mortem data (first marker)"))
NULL