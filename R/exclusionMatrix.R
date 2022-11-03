#' Find the number of incompatible markers for each
#'
#' This function computes the number of exclusions, i.e., the number of
#' incompatible markers, for each pairwise comparison. By default, mutation
#' models are ignored. The main work is done by [forrel::findExclusions()].
#' 
#' @param dvi A `dviData` object, typically created with [dviData()].
#' @param removeMut A logical. If TRUE (default), all mutations models are
#'   stripped.
#'
#' @return An integer matrix with `length(pm)` columns and `length(am)` rows.
#'
#' @examples
#'
#' # Plane crash example
#' exclusionMatrix(planecrash)
#'
#' # Inspect a particular pair: M3 vs V6
#' pm = planecrash$pm
#' am = planecrash$am
#' forrel::findExclusions(am, id = "M3", candidate = pm$V6)
#'
#' # Plot one of the incompatible markers
#' plotPedList(c(am[3], pm[6]), marker = "D7S820", col = list(red = "M3"))
#'
#'
#' @importFrom forrel findExclusions
#' @export
exclusionMatrix = function(dvi, removeMut = TRUE) {
  
  pm = dvi$pm
  am = dvi$am
  missing = dvi$missing
  
  # Ensure that input objects are lists
  if(is.singleton(pm)) 
    pm = list(pm)
  if(is.ped(am)) 
    am = list(am)
  
  npm = length(pm)
  nam = length(am)
  nmiss = length(missing)
  
  # Initialise matrix
  mat = matrix(0L, nrow = npm, ncol = nmiss, dimnames = list(names(pm), missing))
  
  # AM components
  comp = getComponent(am, missing, checkUnique = TRUE, errorIfUnknown = TRUE)
  
  # Loop through each pair of victim vs missing
  for(i in 1:npm) for(j in 1:nmiss) {
    vic = pm[[i]]
    ref = am[[comp[j]]]
    mat[i, j] = length(findExclusions(ref, id = missing[j], candidate = vic))
  }
  
  mat
}
