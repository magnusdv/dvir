
#' @importFrom pedprobr likelihood
loglikTotal = function(x, markers = seq_len(nMarkers(x))) {
  sum(likelihood(x, marker = markers, logbase = exp(1)))
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

ftime = function(st, digits = 3) {
  format(Sys.time() - st, digits = digits)
}

# Fast intersection. NB: assumes no duplicates!
.myintersect = function (x, y) {
  y[match(x, y, 0L)]
}

.mysetdiff = function(x, y) {
  unique.default(x[match(x, y, 0L) == 0L])
}

pluralise = function(noun, n, numberFirst = TRUE) {
  s = if(n == 1) noun else sprintf("%ss", noun)
  if(numberFirst) 
    s = paste(n, s)
  s
}

.safeMax = function(v) {
  if(length(v) == 0 || all(is.na(v))) 
    return(NA)
  max(v, na.rm = TRUE)
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

# Undisputed entries in a LR/GLR matrix
undisputedEntries = function(M, threshold, strict = TRUE) {
  # Replace NAs with -1 to allow numeric comparisons below
  M[is.na(M)] = -1
  
  # Indices of matches exceeding threshold
  highIdx = which(M >= threshold - sqrt(.Machine$double.eps), arr.ind = TRUE)
  
  # Return if empty
  if(!nrow(highIdx)) 
    return(highIdx)
  
  # Identify which rows of `highIdx` to keep
  if(strict) { # undisputed = no others in same row or column exceed 1
    goodRows = which(rowSums(M <= 1) == ncol(M) - 1)
    goodCols = which(colSums(M <= 1) == nrow(M) - 1)
    isUndisp = highIdx[, "row"] %in% goodRows & highIdx[, "col"] %in% goodCols
  }
  else { # undisputed = no others in same row or column exceed LR/threshold
    isUndisp = sapply(seq_len(nrow(highIdx)), function(k) {  # safer than apply(.., 1)!
      rw = highIdx[k,1]
      cl = highIdx[k,2]
      all(c(M[rw, -cl], M[-rw, cl]) <= M[rw, cl]/threshold)
    })
  }

  # Return matrix of indices of undisputed matches   
  highIdx[isUndisp, , drop = FALSE]
}
