#' Add victim samples to a DVI dataset
#'
#' This is a convenience function for adding new victim samples to an existing
#' [dviData()] object. Users are recommended to use this instead of manipulating
#' the `pm` data directly, as it ensures consistency of the output and also
#' recomputes the internal `pairings` object.
#'
#' @param dvi A [dviData()] object.
#' @param newvics A `singleton` object or a list of such.
#'
#' @returns A new [dviData()] object with the new victims added to the `pm`
#'   data.
#'
#' @examples
#' addVictims(example1, singleton("V4"))
#'
#' @export
addVictims = function(dvi, newvics) {
  
  if(is.singleton(newvics))
    newvics = list(newvics)
  else if(!is.pedList(newvics) && pedsize(newvics) == length(newvics))
    stop2("`newvics` must be a `singleton` object or a list of such")
  
  labs = unlist(lapply(newvics, \(x) x$ID), use.names = FALSE)
  if(any(labs %in% names(dvi$pm)))
    stop2("Victim already exists: ", .myintersect(labs, names(dvi$pm)))
  
  # Extend PM data
  newPM = c(dvi$pm, newvics)
  dviData(pm = newPM, am = dvi$am, missing = dvi$missing)
}