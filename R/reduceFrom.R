#' Reduce list of from (PM) samples
#' 
#' @param from A list of singeltons 
#' @param ids.from Character vector with names of PM-samples
#' @export
#' @examples 
#' library(forrel)
#' data(dvi.nfi.mut)
#' dvi.nfi = dvi.nfi.mut
#' from = dvi.nfi$pm
#' from = from[c(1:5, c(3,3,4,5))] #copyies
#' from[[6]]$ID = "V3.c1"
#' from[[7]]$ID = "V3.c2"
#' from[[8]]$ID = "V4.c"
#' from[[9]]$ID = "V5.c"
#' ids.from = c(dvi.nfi$vict, "V3.c1", "V3.c2", "V4.c", "V5.c")
#' red = reduceFrom(from,ids.from)
#' Merge 1, paste navn, fjerne
#' from = from[-(1:length(from))[ids.from == red[1,2]]]
reduceFrom = function(from, ids.from, thresholds = c(MZ = 0.99, PO = 0.99)){
# Estimate kappas for all pairs
  D = t(combn(ids.from, 2))
  n1 = dim(D)[1]
  kappas = IBDestimate(from, D)
# Find candidates for MZ-removal
  linesMZ = (1:n1)[kappas$k2 > thresholds["MZ"]]
  kappas12 = kappas[,1:2]
  remove = NULL
  if(length(linesMZ) > 0) {
  ids = getIDs(from) # in case sorting is not the same as in ids.from
    for (i in 2:length(from)){
      remove = c(remove, (1:(i-1))[ids[i] %in% kappas12[1:(i-1),]])
    }
  }
  
#  linesPO = (1:n1)[kappas$k1 > thresholds["PO"]] #senere
  kappas[linesMZ,]
}


  