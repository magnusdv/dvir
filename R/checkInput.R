#' Checks input and initial likelihood
#' 
#' Parameters are described as in [divSearch], nothing is done here 
#' beyond checking  and calculating NULL likelihoood.
#' 
#' @param from A list of singeltons 
#' @param to A list of pedigrees
#' @param ids.from Character vector with names of PM-samples if these
#' come from different individuals. If NULL, different individuals are not assumed
#' and names are taken from `from`
#' @param ids.to Character vector with names of missing persons
#' @param limit Double. Threshold for LR
#' @param nbest Integer If not `NULL`, search is restricted to the `nbest`
#' marginal solutions for each victim
#' @param extend Logical. If `TRUE`any number of victims need not be among 
#' the missing persons
#' @param plots Logical. If `TRUE`, plot(s) of the best solution will be made. 
#' A plot window should  be prepared and sized in advance
#' 
#' @return List with two elements. First `error` is `NULL` if no error is detected
#' the second `lik0` is the NULL likelihood.
#' 
#' @details In addition to the obvious, it is e.g. checked that all PM-samples have
#' some marker data and that no MP-s have any marker data 
#' @export
#' @examples 
#' \dontrun{
#' data(dvi.nfi)
#' data(dvi.nfi.mut)
#' dvi.nfi = dvi.nfi.mut
#' from = dvi.nfi$pm
#' to = dvi.nfi$am
#' ids.from = dvi.nfi$vict
#' ids.to = dvi.nfi$miss
#' extend = T; limit = -1; nbest = 1; plots = TRUE
#' foo = checkInput(from, to, ids.from, ids.to, extend = extend, limit = limit, 
#'      nbest = nbest, plots = plots)
#' stopifnot(is.null(foo$error) $ foo$lik0 > 0)
#'      
#' v1 = singleton("V1")
#' m1 = marker(v1, V1=1)
#' v1 = addMarkers(v1, m1)
#' v2 = singleton("V2")
#' m2 = marker(v2, V2 = 2)
#' v2 = addMarkers(v2, m2)
#' from = list(v1, v2)
#' ids.from = getIDs(from)
#' to1 = nuclearPed(1, father = "R1", mother = "MO", child = "MP1")
#' m1 = marker(to1, "R1" = 2)
#' to = addMarkers(to1, m1)
#' ids.to = "MP1"
#' foo = checkInput(from, to, ids.from, ids.to, extend = extend, limit = limit, 
#'      nbest = nbest, plots = plots)
#' }

checkInput = function(from, to, ids.from, ids.to, extend = TRUE, limit = 1, 
                      nbest = 1, plots = FALSE){
# Check from: must be a singleton or a list of singeltons.
  if(!is.singleton(from) & !is.pedList(from))
    stop("First argument must be a singleton or a list of such.")
  
# from elements must have some marker data
  if(is.singleton(from) & is.null(from$markerdata))
    stop("No marker data for a PM sample.")
  if(is.pedList(from)){
    noMarkerData =  unlist(lapply(from, function(x) is.null(x$markerdata)))
    if(any(noMarkerData))
      stop("No markerdata for a PM sample.")
  }
  
# Check to
  if(!is.pedList(to) & !is.ped(to))
    stop("Second argument `to`` must be a  ped object or a list of such.")
  if(is.ped(to)){
    if(is.null(to$markerdata))
      stop("All reference families must have at least one genotyped")
  }
  
  # Check ids.from  
  nPM = length(ids.from)
  if(nPM < 1) stop("Names of PM samples must be a vector of length 1 or more.")
  if(length(unique(ids.from)) < nPM )
    errortext = stop("PM samples must be unique.")
  if(!all(sort(ids.from) == sort(getIDs(from))))
      stop("Names of PM samples (ids.from) samples must equal those in from")

# Check ids.to
  nAM = length(ids.to)
  if(nAM < 1) stop("Names of AM samples must be a vector of length 1 or more.")
  if(length(unique(ids.to)) < nAM)
    errortext = stop("AM samples must be unique.")
  if(!all(ids.to %in% getIDs(to)))
    stop("Names of AM samples must be in to")

# Check remaining parameters
  if(!is.logical(extend))
    stop("Parameter extend must be logical")
  if(!is.double(limit))
    stop("Parameter limit must be a double")
  if(nbest < 1) 
    stop("nbest must be NULL or an integer greater than 0")

# Check that all PM-s have marker data and that no MP-s have
    g = getAlleles(from)
    allTRUE = apply(g, 1, function(x) any(!is.na(x)))
    if(!all(allTRUE))
      stop(paste("Some PM-s without marker data:", ids.from[!allTRUE]))
    g = getAlleles(to, ids = ids.to)
    if(!all(is.na(g)))
      stop("Genotypes for missing persons.")

# Calculate initial likelihood    
    lik0 = prod(LR(list(from, to), 1)$likelihoodsPerSystem)
    if (lik0 == 0){
      stop("Initial data has 0 likelihood.")
   }

  list(error = NULL, lik0 = lik0)
}