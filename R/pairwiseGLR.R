pairwiseGLR = function(dvi, jointTable = NULL, LRmatrix = NULL, threshold = 1e4, verbose = FALSE) {
  
  dvi = consolidateDVI(dvi)
  pm = dvi$pm
  am = dvi$am
  vics = names(pm)
  missing = dvi$missing
  
  if(threshold <= 1) 
    stop2("Argument `threshold` must exceed 1: ", threshold)
  
  # Pairwise LR matrix (as additional info)
  LRmat = LRmatrix %||% pairwiseLR(dvi, verbose = FALSE)$LRmatrix
  
  # Joint table (should normally be supplied)
  joint = jointTable %||% dviJoint(dvi, verbose = verbose)
  
  # Initialise GLR matrix
  G = matrix(NA_real_, nrow = length(vics), ncol = length(missing), 
             dimnames = list(vics, missing))
  
  # Fill matrix row-wise
  for (v in vics) {
    # corresponding column of joint table
    a = joint[[v]]
    
    # If all equal (to '*' or some Mj): GLR not defined
    if(length(unique.default(a)) < 2)
      next
    
    # Index of 1st row with v = m, for each m in miss
    idxNumer = match(missing, a)
    
    # Index of 1st row with v != m
    idxDenom = ifelse(idxNumer == 1, which(a != a[1])[1], 1)
    
    # Insert in matrix
    ll = joint$loglik
    G[v,  ] = exp(ll[idxNumer] - ll[idxDenom])
  }
  
  # Indices of matches exceeding GLR threshold
  highGLR = which(G >= threshold - sqrt(.Machine$double.eps), arr.ind = TRUE)
  
  foundVics = vics[highGLR[, 1]]
  foundMiss = missing[highGLR[, 2]]
  
  # Summary
  summ = NULL
  
  if(length(foundVics)) {
    summ = data.frame(
      Family = getFamily(dvi, foundMiss),
      Missing = foundMiss, 
      Sample = foundVics, 
      LR = LRmat[highGLR],
      GLR = G[highGLR],
      Conclusion = "Match (GLR)",
      Comment = sprintf("Joint analysis {%s}", paste0(missing, collapse = ",")),
      row.names = NULL)
  }
  
  # Data from identified samples (keep for updating dviRed below)
  foundData = pm[foundVics]
   
  # Reduced DVI dataset
  newvics = setdiff(vics, foundVics)
  newmissing = setdiff(missing, foundMiss)
  dviRed = subsetDVI(dvi, pm = newvics, missing = newmissing, verbose = verbose)
  
  # Move vic data to AM data
  names(foundVics) = foundMiss
  relevantMP = intersect(foundMiss, labels(am))
  if(length(relevantMP))
    dvi$am = transferMarkers(from = foundData, to = am, 
                             idsFrom = foundVics[relevantMP], 
                             idsTo = relevantMP, erase = FALSE)
    
  list(dviReduced = dviRed, GLRmatrix = G, summary = summ)
}