#' Combinations without duplications
#'
#' This is similar to [expand.grid()] except that combinations with repeated elements are not
#' included. The element "*" is treated separately, and is allowed to be repeated.
#'
#' `gridSize()` calculates the number of rows produced by [expand.grid.nodup()], either exactly or
#' as a guaranteed lower bound. Its main purpose is to detect and abort cases that are too large.
#' The exact algorithm counts matchings in the corresponding bipartite row-value graph.
#'
#' @param lst A list of vectors.
#' @param max A positive integer. If the number of combinations exceeds this, the function aborts
#'   with an informative error message. Default: 1e5.
#' @param verbose A logical.
#'
#' @return A data frame.
#' @seealso [expand.grid()], [ncomb()]
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
expand.grid.nodup = function(lst, max = 1e5, verbose = FALSE) {
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
  
  # Exit with error if too many
  b = gridSize(lst, max = max)
  if(verbose)
    cat(sprintf("Grid size: %d (%s)\n", b$size, b$type))
  if(b$size > max)
    stop2(sprintf("Grid size (%s: %g) exceeds max (%d)", b$type, b$size, max))
    
  #msg2 = "Possible strategies:\n * Increase `limit`\n * Decrease `threshold`\n * Increase `maxAssign`\n     #* Set `strict` to FALSE"
  #stop2(paste0(msg1, msg2))
    
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

#' @rdname expand.grid.nodup
gridSize = function(lst, max = Inf) {
  # first quick exact count, then exact if not too big, otherwise lower bound only.
  b = gridSizeExact(lst, cap = 10001)
  if(b <= 10000) 
    return(list(size = b, type = "exact"))
  
  b = gridSizeLowerBound(lst)
  if(b > max) 
    return(list(size = b, type = "lower bound"))
  
  b = gridSizeExact(lst, cap = max + 1)
  if(b < max) 
    return(list(size = b, type = "exact"))
  else 
    return(list(size = b, type = "lower bound"))
}
  
gridSizeExact = function(lst, cap = Inf) {
  # Exact count, returned as min(count, cap).
  
  n = length(lst)
  len = lengths(lst)
  
  if(n == 0L || any(len == 0L))
    return(0)
  
  row = rep.int(seq_len(n), len)
  val = as.character(unlist(lst, use.names = FALSE))
  star = val == "*"
  
  local = tabulate(row[star], nbins = n)
  row = row[!star]
  val = val[!star]
  
  if(!length(val))
    return(min(prod(local), cap))
  
  # Encode ordinary values as integers
  col = match(val, unique.default(val))
  
  # Values occurring in one row only are private (and hence local)
  private = tabulate(col, nbins = max(col))[col] == 1L
  
  if(any(private)) {
    local = local + tabulate(row[private], nbins = n)
    row = row[!private]
    col = col[!private]
  }
  
  if(!length(col))
    return(min(prod(local), cap))
  
  # Re-index shared values to 1,...,m
  col = match(col, unique.default(col))
  m = max(col)
  
  if(m > 53L)
    stop2("Too many shared values; exact counting not supported")
  
  # Row-indexed list of remaining shared values.
  vals = split(col, factor(row, levels = seq_len(n)))
  deg = lengths(vals)
  
  # If any row has no choices at all, the grid is empty
  if(any(deg + local == 0L))
    return(0)
  
  # Constrained rows first; rows with local choices are less constrained
  ord = order(local > 0L, deg + local)
  local = local[ord]
  vals = vals[ord]
  
  bit = 2^(seq_len(m) - 1L)
  bits = lapply(vals, \(v) bit[v])
  
  # Memoization environments for rows 1,...,n; keys are bit-masks of used shared values
  memo = replicate(n + 1L, new.env(parent = emptyenv(), hash = TRUE), simplify = FALSE)
  
  # Recursively count completions of rows i,...,n 
  rec = function(i, mask) {
    if(i > n)
      return(1)
    
    key = as.character(mask)
    old = memo[[i]][[key]]
    if(!is.null(old))
      return(old)
    
    # Local choice: "*" or a private ordinary value.
    total = if(local[i]) local[i]*rec(i + 1L, mask) else 0
    
    # Shared choice: only allowed if the corresponding bit is unused.
    if(total < cap) {
      bi = bits[[i]]
      
      for(b in bi[floor(mask/bi) %% 2 == 0]) {
        total = total + rec(i + 1L, mask + b)
        
        if(total >= cap) {
          total = cap
          break
        }
      }
    }
    
    memo[[i]][[key]] = total
    total
  }
  
  rec(1L, 0)
}


gridSizeLowerBound = function(lst, cap = Inf) {
  n = length(lst)
  len = lengths(lst)
  
  if(n == 0L || any(len == 0L))
    return(0)
  
  row = rep.int(seq_len(n), len)
  val = as.character(unlist(lst, use.names = FALSE))
  isStar = val == "*"
  
  # local[i] = choices in row i that cannot block anything: "*" + private ordinary values
  local = tabulate(row[isStar], nbins = n)
  row = row[!isStar]
  val = val[!isStar]
  
  if(!length(val))
    return(min(prod(local), cap))
  
  # Encode ordinary values as integers; no within-row duplicates assumed
  col = match(val, unique.default(val))
  
  # Ordinary values occurring in only one row are private, hence local
  private = tabulate(col, nbins = max(col))[col] == 1L
  
  if(any(private)) {
    local = local + tabulate(row[private], nbins = n)
    row = row[!private]
    col = col[!private]
  }
  
  if(!length(row))
    return(min(prod(local), cap))
  
  vals = split.default(col, factor(row, levels = seq_len(n)))
  deg = lengths(vals)
  
  if(any(local == 0L & deg == 0L))
    return(0)
  
  # Empirically best cheap order: few total choices first
  ord = order(deg + local, method = "radix")
  seen = rep(FALSE, max(col))
  dp = 1
  
  for(i in ord) {
    active = which(dp > 0)
    if(!length(active))
      return(0)
    
    t = active - 1L
    ndp = numeric(length(dp) + 1L)
    
    # Local choices: t unchanged
    if(local[i] > 0L)
      ndp[active] = dp[active] * local[i]
    
    # Shared choices: t increases by 1; only previously seen values can be blocked
    vi = vals[[i]]
    
    if(length(vi)) {
      overlap = sum(seen[vi])
      avail = deg[i] - pmin.int(t, overlap)
      ok = avail > 0
      if(any(ok))
        ndp[active[ok] + 1L] = ndp[active[ok] + 1L] + dp[active[ok]]*avail[ok]
      
      seen[vi] = TRUE
    }
    
    # Cap to avoid overflow; does not affect the final bound.
    ndp[ndp > cap] = cap
    
    dp = ndp
  }
  
  min(sum(dp), cap)
}

# For testing: random pairing list of length n, with m vics
randomPairingList = function(n, m) {
  vals = paste0("V", seq_len(m))
  stars = as.logical(sample.int(2L, n, replace = TRUE) - 1)
  lens = sample.int(m, n, replace = TRUE)
  lapply(1:n, \(i) c(if(stars[i]) "*", sample(vals,lens[i])))
}

