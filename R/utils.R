
#' @importFrom pedprobr likelihood
loglikTotal = function(x, markers = seq_len(nMarkers(x))) {
  sum(likelihood(x, marker = markers, logbase = exp(1), eliminate = 1))
}


# Modified version of stop()
stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}


`%||%` = function(x, y) {
  if(is.null(x)) y else x
}

trunc = function(x, printMax = 10) {
  if(length(x) <= printMax)
    return(toString(x))
  y = c(x[1:5], "...", x[length(x)])
  toString(y)
}

# example: sprintfNamed("I am %{name}s", name = "your father")
sprintfNamed = function(fmt, ...) {
  arglist = list(...)
  if(length(arglist) == 1 && is.list(arglist[[1]]))
    arglist = arglist[[1]]
  
  nms = names(arglist)
  if(is.null(nms) || any(nms == "")) 
    stop2("Arguments must have names")
  
  patterns = sprintf("%%{%s}", nms)
  hasPat = vapply(patterns, function(p) grepl(p, fmt, fixed = TRUE), FALSE)
  
  # Remove unused arguments
  arglist = arglist[hasPat]
  patterns = patterns[hasPat]

  for(i in seq_along(patterns)) {
    fmt = gsub(pattern = patterns[i],
               replacement = sprintf("%%%d$", i),
               fmt, fixed = TRUE)
  }

  miss = regmatches(fmt, m = gregexec("%\\{([^%]*)}", fmt, perl = TRUE))[[1]]
  if(length(miss))
    stop2("Missing values for variable: ", sprintf("'%s'", miss[2, ]))
  
  do.call(sprintf, append(arglist, fmt, 0))
}

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