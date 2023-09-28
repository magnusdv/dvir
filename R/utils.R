
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
