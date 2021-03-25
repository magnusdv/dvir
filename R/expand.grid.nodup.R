#' Combinations without duplications
#'
#' This is a simple extension of [expand.grid()] which removes all combinations
#' with repeated elements.
#'
#' @param lst A list of vectors.
#'
#' @author Magnus Dehli Vigeland
#'
#' @return A data frame.
#' @seealso [expand.grid()]
#'
#' @examples
#'
#' lst = list(1:2, 1:2)
#'
#' # Compare
#' expand.grid.nodup(lst)
#' expand.grid(lst)
#'
#'
#' @export
expand.grid.nodup = function(lst) {
  if(is.data.frame(lst))
    stop("Unexpected input: Argument `lst` is a data frame")
  if(!is.list(lst))
    stop("Argument `lst` should be a list")
  
  len = length(lst)  
  if(len == 0)
    return(data.frame())
  
  nms = names(lst)
  if(is.null(nms)) 
    nms = paste0("Var", 1:len)
  
  # If empty: Return early
  if(any(lengths(lst) == 0))
    return(data.frame(matrix(nrow = 0, ncol = len), dimnames = list(NULL, nms)))
  
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
  
  # Convert to data frame and add names
  res = as.data.frame.matrix(do.call(rbind, reslist))
  names(res) = nms
  
  res
}

