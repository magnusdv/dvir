#' Combinations without duplications
#'
#' This is a simple extension of [expand.grid()] which removes all combinations with repeated elements.
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
  
  # Grid
  args = c(lst, list(KEEP.OUT.ATTRS = F, stringsAsFactors = F))
  full.grid = do.call(expand.grid, args)
  
  # Remove rows with repeated elements
  nodups = apply(full.grid, 1, anyDuplicated.default, incomparables = "*") == 0
  res = full.grid[nodups, , drop = F]
  
  # Fix row names
  rownames(res) = NULL
  
  res
}
