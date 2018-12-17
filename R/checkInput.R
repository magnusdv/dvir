#' @export
checkInput = function(from, to, ids.from , ids.to ){
  errortext = NULL
  if(!pedtools::is.pedList(from) & !pedtools::is.ped(from))
    errortext = c("First argument must be a  ped object or a list of such")
  if(!pedtools::is.pedList(to) & !pedtools::is.ped(to))
    errortext = c("Second argument must be a  ped object or a list of such")
  nvp = length(ids.from)
  if(nvp < 1) errortext = c("Third argument must be a vector of length 1 or more")
  nmp = length(ids.to)
  if(nmp < 1) errortext = c("Four argument must be a vector of length 1 or more")
  if(pedtools::is.ped(from)){
    labelspm = labels(from)
    if(!all(ids.from %in% labels(from)))
      errortext = c("Third argument must be a subset of labels of first argument")
  }
  if(pedtools::is.ped(to)){
    labelsam = labels(to)
    if(!all(ids.to %in% labels(to)))
      errortext = c("Fourth argument must be a subset of labels of second argument")
  }
  if(pedtools::is.pedList(from)){
    labelspm = unlist(lapply(from, labels))
    if(!all(ids.from %in% labelspm))
      errortext = c("Third argument must be a subset of labels of first argument")
  }
  if(pedtools::is.pedList(to)){
    labelsam = unlist(lapply(to, labels))
    if(!all(ids.to %in% labelsam))
      errortext = c("Fourth argument must be a subset of labels of second argument")
  }
  if(length(unique(ids.from)) < nvp )
    errortext = c("Victims must be unique")
  if(length(unique(ids.to)) < nmp )
    errortext = c("Missing persons must be unique")
  g = getAlleles(from)
  if(all(is.na(g)))
    errortext = c("No genotypes for victims")
  g = getAlleles(to, ids = ids.to)
  if(!all(is.na(g)))
    errortext = c("Genotypes for missing persons")
  likNULL = prod(LR(list(from, to), 1)$likelihoodsPerSystem)
  if (likNULL == 0){
    errortext = "Initial data has 0 likelihood"
   }
  list(error = errortext, lik0 = likNULL)
}