#' Extract the marker database from a DVI dataset
#'
#' @param dvi A [dviData()] object.
#' @param check A logical. If TRUE (default), the functions checks if the AM and
#'   PM databases are identical, and raises an error if not. If FALSE, the AM
#'   database is returned without checking PM.
#' @param includeAttrs A logical, indicating if the full list of attributes of
#'   each marker (including mutation model and position data) should be
#'   returned, or only the allele frequencies. Default: FALSE.
#'
#' @returns A list of marker information.
#'
#' @examples
#' getDatabase(example1)
#' getDatabase(example1, includeAttrs = TRUE)
#'
#' @export
getDatabase = function(dvi, check = TRUE, includeAttrs = FALSE) {
  dvi = consolidateDVI(dvi)
  getFUN = if(includeAttrs) getLocusAttributes else getFreqDatabase
  
  dbAM = getFUN(dvi$am) 
  if(!check) 
    return(dbAM)
  
  dbPM = getFUN(dvi$pm)
  
  if(!identical(dbAM, dbPM))
    stop2("Databases differ between AM and PM")
  
  dbAM
}
