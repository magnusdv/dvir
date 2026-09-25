#' Set or add victim samples
#'
#' Convenience functions for replacing or adding victim samples in an existing [dviData()]
#' object. These functions should generally be preferred over direct manipulation of the
#' `pm` component, as they ensure consistency and recompute the internal `pairings`
#' object.
#'
#' `setVictims()` replaces all existing victim samples, while `addVictims()` retains the
#' existing samples and appends new ones.
#'
#' @param dvi A [dviData()] object.
#' @param vics A `singleton` object or a list of such.
#'
#' @returns A modified [dviData()] object.
#'
#' @examples
#' # Add a victim
#' addVictims(example1, singleton("V4"))
#' 
#' # Replace the PM data
#' setVictims(example1, singleton("V4"))
#'
#' @export
setVictims = function(dvi, vics) {
  if(is.singleton(vics))
    vics = list(vics)
  else if(!all(vapply(vics, is.singleton, logical(1))))
    stop2("`vics` must be a `singleton` object or a list of such")

  dviData(pm = vics, am = dvi$am, missing = dvi$missing)
}

#' @rdname setVictims
#' @export
addVictims = function(dvi, vics) {
  
  if(is.singleton(vics))
    vics = list(vics)
  else if(!all(vapply(vics, is.singleton, logical(1))))
    stop2("`vics` must be a `singleton` object or a list of such")
  
  labs = unlist(lapply(vics, \(x) x$ID), use.names = FALSE)
  if(any(labs %in% names(dvi$pm)))
    stop2("Victim already exists: ", .myintersect(labs, names(dvi$pm)))
  
  # Extend PM data
  newPM = c(dvi$pm, vics)
  dviData(pm = newPM, am = dvi$am, missing = dvi$missing)
}
