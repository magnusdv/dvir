#' Export DVI data to Familias
#'
#' Convenience wrapper for [pedFamilias::writeFam()] with `dvi` parameter set to TRUE.
#'
#' @param x A `dviData` object.
#' @inheritParams pedFamilias::writeFam
#'
#' @return The file name is returned invisibly.
#'
#' @seealso [familias2dvir()], [pedFamilias::writeFam()]
#'
#' @examples
#' writeFamDVI(planecrash, tempfile(fileext = ".fam"))
#'
#' @export
writeFamDVI = function(x, famfile = "dvi.fam", params = NULL, dbOnly = FALSE,
                        openFam = FALSE, FamiliasPath = NULL, verbose = TRUE) {
  x = consolidateDVI(x)
  
  # Ensure DVI mode
  params$dvi = TRUE
  
  pedFamilias::writeFam(x, famfile = famfile, params = params, dbOnly = dbOnly,
                        openFam = openFam, FamiliasPath = FamiliasPath, 
                        verbose = verbose)
}
