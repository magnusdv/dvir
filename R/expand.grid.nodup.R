#' Combinations without duplications
#'
#' This is similar to [expand.grid()] except that combinations with repeated
#' elements are not included. The element "*" is treated separately, and is
#' allowed to be repeated.
#'
#' @param lst A list of vectors.
#' @param max A positive integer. If the number of combinations (according to a
#'   preliminary lower bound) exceeds this, the function aborts with an
#'   informative error message. Default: 1e5.
#'
#' @return A data frame.
#' @seealso [expand.grid()]
#'
#' @examples
#'
#' lst = list(1, 1:2, 3:4)
#'
#' # Compare
#' expand.grid.nodup(lst)
#' expand.grid(lst)
#'
#' # Typical use case for DVI
#' lst2 = generatePairings(example1)
#' expand.grid.nodup(lst2)
#'
#' @export
expand.grid.nodup = function(lst, max = 1e5) {
  if(is.data.frame(lst))
    stop2("Unexpected input: Argument `lst` is a data frame")
  if(!is.list(lst))
    stop2("Argument `lst` should be a list")
  
  len = length(lst)  
  if(len == 0)
    return(data.frame())
  
  nms = names(lst)
  if(is.null(nms)) 
    nms = paste0("Var", 1:len)
  
  # If empty: Return early
  if(any(lengths(lst) == 0))
    return(data.frame(matrix(nrow = 0, ncol = len, dimnames = list(NULL, nms))))
  
  
  # Warn if too many
  b = estimateGridSize(lst)
  if(b > max) {
    msg1 = sprintf("\nNumber of assignments is at least %g, exceeding the current limit of %d.\n", b, max)
    msg2 = "Possible strategies:\n\n * Increase `limit`\n * Decrease `threshold`\n * Increase `maxAssign`\n * Set `strict` to FALSE"
    stop2(paste0(msg1, msg2))
  }
  
  # Recursive function
  recurse = function(a) {
    if(length(a) == 1)
      return(as.list.default(a[[1]]))
    
    r = lapply(a[[1]], function(val) {
      b = a[-1]
      if(val != "*")
        b = lapply(b, function(vec) vec[vec != val])
      lapply(recurse(b), function(vec) c(val, vec))
    })
    
    unlist(r, recursive = FALSE)
  }
  
  # Call it!
  reslist = recurse(lst)
  if(!length(reslist))
    stop2("No possible solutions (empty grid)")
  
  # Convert to data frame and add names
  resmat = matrix(unlist(reslist, recursive = FALSE, use.names = FALSE), 
                  byrow = TRUE, ncol = len)
  res = as.data.frame.matrix(resmat)
  names(res) = nms
  
  res
}


# Compute a lower bound on the number of assignments for a given list of pairings
# Exact if entries of lst are all contained in each other
estimateGridSize = function(lst) {
  
  n = length(lst)
  lens = sort.int(lengths(lst))
  
  if(n == 0 || any(lens == 0))
    return(0)
  
  allstars = all(vapply(lst, function(v) length(v) > 0 && "*" %in% v, TRUE))
  
  if(!allstars)
    # Assume no stars. Trivial formula (Theorem 2 from paper 2022)
    return(prod(lens - 0:(n-1)))
  
  # From here: All sets contain "*"
  C = lens - 1
  a = matrix(0, ncol = n + 1, nrow = n)
  
  # first column: a_n,0 (no stars)
  a[, 1] = cumprod(C - 0:(n-1))
  
  # Diagonal: k stars
  for(k in 1:n) 
    a[k, k+1] = 1
  
  # Remaining entries: Recurse!
  for(k in seq_len(n)) for(s in seq_len(k-1))
    a[k,s+1] = a[k-1, s] + a[k-1, s+1] * (C[k] - (k-1-s))
  
  sum(a[n, ])
}
