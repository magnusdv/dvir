#' Combine summary tables
#'
#' Combines summary tables from various functions into a final result table.
#'
#' The output is controlled by `centricity`, which the following effect:
#' 
#' * `AM`: 
#'     - Column order: `Family`, `Missing`, `Sample`, `LR`, `GLR`, `Conclusion`, `Comment`
#'     - Sort by: `Family` and `Missing`
#' * `PM`:
#'     - Column order: `Sample`, `Missing`, `Family`, `LR`, `GLR`, `Conclusion`, `Comment`
#'     - Sort by: `Sample`
#' 
#' Columns (in any of the data frames) other than these are simply ignored.
#' 
#' @param dfs A list of data frames.
#' @param centricity Either "AM" or "PM", controlling column order and sorting.
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
#' combineSummaries(list(u$summary, a$summary$AM))
#' combineSummaries(list(u$summary, a$summary$PM), centricity = "PM", dvi = planecrash)
#' 
#' @export
combineSummaries = function(dfs, centricity = c("AM", "PM"), dvi = NULL) {
  centr = match.arg(centricity)
  orderBy = switch(centr, AM = c("Family", "Missing"), PM = "Sample")
  
  # Column order depends on centricity
  allCols = switch(centr,
    AM = c("Family", "Missing", "Sample", "LR", "GLR", "Conclusion", "Comment"),
    PM = c("Sample", "Missing", "Family", "LR", "GLR", "Conclusion", "Comment")
  )
  
  # Harmonize data frames to have all columns
  dfsExt = lapply(dfs, function(df) {
    if(is.null(df) || nrow(df) == 0) 
      return(NULL)
    
    # Missing columns
    miscol = !allCols %in% names(df)
    if(all(miscol))
      return(NULL)
    
    # Add missing columns; fill with NAs
    for(cc in allCols[miscol])
      df[[cc]] = NA
    
    df[allCols]
  })
  
  final = do.call(rbind, dfsExt)
  
  # If no ordering: Return
  if(is.null(orderBy))
    return(final)
  
  # If no `dvi` object provided, order in standard way
  if(is.null(dvi)) {
    final = final[do.call(order, final[orderBy]), , drop = FALSE]
    return(final)
  }
  
  # Take ordering data from DVI object
  ordvec = lapply(orderBy, function(cc) {
    switch(cc,
     Family = match(final$Family, names(dvi$am)),
     Missing = match(final$Missing, dvi$missing),
     Sample = match(final$Sample, names(dvi$pm)),
     stop2("Column does not exist: ", cc))
  })
  
  final = final[do.call(order, ordvec), , drop = FALSE]
  rownames(final) = NULL
  final
}