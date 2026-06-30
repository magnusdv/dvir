#' Assignment grid for joint analysis
#'
#' Given a list of candidate pairings for each victim, `expand.grid.nodup()`
#' constructs the corresponding assignment table. This is similar to [expand.grid()],
#' but excludes assignments where the same missing person is used more than once. 
#' The special value "*" is treated as non-assignment and may be repeated. 
#' This is used by [dviJoint()] to enumerate joint assignments.
#'
#' @param pairings A list of character vectors, typically produced by
#'   [generatePairings()]. Each list entry gives the possible assignments for
#'   one victim.
#' @param max A positive number. If the grid size exceeds this, the function
#'   aborts with an informative error message. Default: 1e5.
#' @param df A logical, indicating if the function should return a data frame
#'   (default) or a matrix.
#' @param verbose A logical.
#' @param lst Deprecated alias for `pairings`.
#'
#' @return A data frame or matrix with one row for each valid assignment.
#'
#' @seealso [dviGridSize()], [dviJoint()], [generatePairings()], [expand.grid()]
#'
#' @examples
#' pairings = list(V1 = c("M1", "M2", "*"), V2 = c("M1", "M2", "*"))
#' expand.grid.nodup(pairings)
#' 
#' # Compare with expand.grid()
#' expand.grid(pairings)
#' 
#' # Dataset `example1`
#' pairings = generatePairings(example1)
#' g = expand.grid.nodup(pairings)
#' n = dviGridSize(pairings = pairings)
#' stopifnot(nrow(g) == n$size)
#'
#' @export
expand.grid.nodup = function(pairings = NULL, max = 1e5, df = TRUE, verbose = FALSE, lst = NULL) {
  if(!is.null(lst)) {
    if(!is.null(pairings))
      stop2("Use only one of `pairings` and `lst`")
    pairings = lst
  }

  lst = pairings

  if(is.data.frame(lst))
    stop2("Unexpected input: Argument `lst` is a data frame")
  if(!is.list(lst))
    stop2("Argument `lst` should be a list")
  
  len = length(lst)  
  if(len == 0)
    return(if(df) data.frame() else matrix(ncol = 0, nrow = 0))
  
  nms = names(lst)
  if(is.null(nms)) 
    nms = paste0("Var", 1:len)
  
  # If empty: Return early
  if(any(lengths(lst) == 0)) {
    emptymat = matrix(nrow = 0, ncol = len, dimnames = list(NULL, nms))
    return(if(df) as.data.frame.matrix(emptymat) else emptymat)
  }
   
  # Exit with error if too many
  b = dviGridSize(pairings = lst, max = max)
  if(verbose)
    cat(sprintf("Grid size: %d (%s)\n", b$size, b$type))
  if(b$size > max)
    stop2(sprintf("Grid size (%s: %g) exceeds max (%d)", b$type, b$size, max))
    
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
  res = matrix(unlist(reslist, recursive = FALSE, use.names = FALSE), 
               byrow = TRUE, ncol = len)
  colnames(res) = nms
  if(df)
    res = as.data.frame.matrix(res)
  
  res
}



