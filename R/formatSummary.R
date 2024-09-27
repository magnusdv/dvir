#' Format final summary table
#'
#' Combines and harmonises summary tables from different DVI analyses
#'
#' The default column order is controlled by `orientation`, which the following
#' effect:
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
#' @param orientation Either "AM" or "PM", controlling column order and sorting.
#' @param columns A (optional) character vector with column names in the wanted
#'   order.
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
#' formatSummary(list(u$summary, a$summary$AM))
#' formatSummary(list(u$summary, a$summary$PM), orientation = "PM", dvi = planecrash)
#'
#' @export
formatSummary = function(dfs, orientation = c("AM", "PM"), columns = NULL, dvi = NULL) {
  
  centr = match.arg(orientation)
  
  # Column order: depends on orientation
  allCols = columns %||% switch(centr,
    AM = c("Family", "Missing", "Sample", "LR", "GLR", "Conclusion", "Comment"),
    PM = c("Sample", "Missing", "Family", "LR", "GLR", "Conclusion", "Comment")
  )
  
  # Harmonize data frames to have all columns
  dfsExt = lapply(dfs, function(df) {
    if(is.null(df) || (nrow(df) == 0)) 
      return(NULL)
    
    # AM must have `Missing`; PM must have `Sample`
    df = switch(centr,
      AM = df[!is.na(df$Missing), , drop = FALSE],
      PM = df[!is.na(df$Sample), , drop = FALSE]
    )
    
    if(nrow(df) == 0) 
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
  
  # If no `dvi` object provided, order in standard way
  if(is.null(dvi)) {
    if(!is.null(final)) {
      orderBy = switch(centr, AM = c("Family", "Missing"), PM = "Sample")
      final = final[do.call(order, final[orderBy]), , drop = FALSE]
    }
    return(final)
  }
  
  fullDf = switch(orientation, 
      AM = data.frame(Family = getFamily(dvi, dvi$missing), Missing = dvi$missing, row.names = NULL),
      PM = data.frame(Sample = names(dvi$pm), row.names = NULL))
  
  # If no data, return fullDf
  if(is.null(final)) {
    for (cc in setdiff(allCols, names(fullDf)))
      fullDf[[cc]] = NA
    return(fullDf)
  }

  orderBy = switch(centr, AM = "Missing", PM = "Sample")
  final = merge(fullDf, final, all.x = TRUE, sort = FALSE)
  final = final[match(fullDf[[orderBy]], final[[orderBy]]), , drop = FALSE]
  rownames(final) = NULL
  final
}
