#' Generate list of assignments
#'
#' Based on a list or a matrix describiing possible moves of victims to missing
#' persons, all unique possibilities are generated.
#' 
#' @param lst List or 0-1 matrix describing possible moves.
#' 
#' @author Magnus Dehli Vigeland
#' 
#' @return List of data frames.
#' @export
#' 
#' @examples
#' 
#' # Example 1
#' moves = list(V1 = c("V1", "MP1"), V2 = c("V2", "MP2"))
#' expand.grid.nodup(moves)
#' expand.grid.nodup(moves[1])
#'
#' # Example 2
#' expand.grid.nodup(list(c(1,2,4), c(1,3), 2))
#' 
#' # Example 3
#' expand.grid.nodup(list(c(12, 1), 2))
#' 
#' \dontrun{
#' lst = rep(list(1:12), 5)
#' gr = expand.grid.nodup(lst)   # takes a few seconds
#' length(gr)   # 95040
#' }
#' 
expand.grid.nodup = function(lst) {
  stopifnot(is.list(lst))

  # Grid
  args = c(lst, list(KEEP.OUT.ATTRS = F, stringsAsFactors = F))
  full.grid = do.call(expand.grid, args)

  # Which rows have no duplications?
  nodups = apply(full.grid, 1, anyDuplicated) == 0
  #if(dim(full.grid)[1] == 1) nodups = TRUE # tried removed Thore July 11
  # Good combinations:
  good.grid = full.grid[nodups, , drop = F]

  # Convert from data frame to list
  nr = nrow(good.grid)
  lapply(seq_len(nr), function(i) good.grid[i, , drop = F])
}
