#' Combine summary tables
#'
#' Combines summary tables from various functions into a final result table.
#'
#' @param dfs A list of data frames.
#' @param orderBy A character with column names to sort by.
#' @param dvi A `dviData` object used for sorting. Note that if given, this must
#'   contain all victims and families.
#'
#' @return A data frame.
#'
#' @examples
#' u = findUndisputed(planecrash)
#' a = amDrivenDVI(u$dviReduced, threshold2 = 500)
#'
#' u$summary
#' a$summary
#'
#' combineSummaries(list(u$summary, a$summary),
#'                  orderBy = c("Family", "Missing"),
#'                  dvi = planecrash)
#' @export
combineSummaries = function(dfs, orderBy = NULL, dvi = NULL) {
  allCols = c("Family", "Missing", "Sample","LR", "Conclusion", "Comment")
  #unique(unlist(lapply(dfs, names)))
  
  # Harmonize data frames to have all columns
  dfsExt = lapply(dfs, function(df) {
    if(is.null(df) || nrow(df) == 0) 
      return(NULL)
    
    # Missing columns
    miscol = setdiff(allCols, names(df))
    for(cc in miscol)
      df[[cc]] = NA
    
    df[allCols]
  })
  
  final = do.call(rbind, dfsExt)
  if(is.null(orderBy))
    return(final)
  
  if(is.null(dvi)) {
    final = final[do.call(order, final[orderBy]), , drop = FALSE]
    return(final)
  }
  
  # Use ordering in provided DVI object
  ordvec = lapply(orderBy, function(cc) switch(cc,
                                               Family = match(final$Family, names(dvi$am)),
                                               Missing = match(final$Missing, dvi$missing),
                                               Sample = match(final$Sample, names(dvi$pm)),
                                               stop2("Column does not exist: ", cc)))
  
  final = final[do.call(order, ordvec), , drop = FALSE]
  rownames(final) = NULL
  final
}