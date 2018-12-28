#' Checks input
#' 
#' @param from A list of singeltons 
#' @param to A list of pedigrees
#' @param ids.from Character vector with names of victims
#' @param ids.to Character vector with names of missing persons
#' @param moves A data frame with columns giving the possible moves
#' @return List with two elements. First is `NULL` if no error is detected
#' and otherwise an error text.The other is the NULL likelihood
#' 
#' @export
#' @examples 
#' checkInput(1, 2, 3, 4)
#' 
checkInput = function(from, to, ids.from , ids.to, moves = NULL ){
  errortext = NULL
  if(!pedtools::is.pedList(from) & !pedtools::is.ped(from))
    errortext = c("First argument must be a  ped object or a list of such.")
  if(!pedtools::is.pedList(to) & !pedtools::is.ped(to))
    errortext = paste( errortext,
                c("Second argument must be a  ped object or a list of such."))
  nvp = length(ids.from)
  if(nvp < 1) errortext = paste(errortext, 
              c(" Third argument must be a vector of length 1 or more."))
  nmp = length(ids.to)
  if(nmp < 1) errortext = paste(errortext, 
              c(" Four argument must be a vector of length 1 or more."))
  if(pedtools::is.ped(from)){
    labelspm = labels(from)
    if(!all(ids.from %in% labels(from)))
      errortext = paste(errortext, 
      c(" Third argument must be a subset of labels of first argument."))
  }
  if(pedtools::is.ped(to)){
    labelsam = labels(to)
    if(!all(ids.to %in% labels(to)))
      errortext = paste(errortext, 
      c(" Fourth argument must be a subset of labels of second argument."))
      }
  if(pedtools::is.pedList(from)){
    labelspm = unlist(lapply(from, labels))
    if(!all(ids.from %in% labelspm))
      errortext = paste(errortext,
      c(" Third argument must be a subset of labels of first argument."))
  }
  if(pedtools::is.pedList(to)){
    labelsam = unlist(lapply(to, labels))
    if(!all(ids.to %in% labelsam))
      errortext = paste(errortext, 
      c(" Fourth argument must be a subset of labels of second argument."))
  }
  if(length(unique(ids.from)) < nvp )
    errortext = paste(errortext, c(" Victims must be unique."))
  if(length(unique(ids.to)) < nmp )
    errortext = paste(errortext, c(" Missing persons must be unique."))
  if(!is.null(moves)){
    if(!is.data.frame(moves))
      errortext = paste(errortext, "Last argument must be a dataframe.")
    if(dim(moves)[2] != 2)
      errortext = paste(errortext, "Last argument must be a dataframe with two columns")
  }
  if(!is.null(errortext))
    stop(errortext)
  g = getAlleles(from)
  if(all(is.na(g)))
    errortext = paste(errortext, c("No genotypes for victims."))
  g = getAlleles(to, ids = ids.to)
  if(!all(is.na(g)))
    errortext = paste(errortext, c("Genotypes for missing persons."))
  lik0 = prod(LR(list(from, to), 1)$likelihoodsPerSystem)
  if (lik0 == 0){
    errortext = paste(errortext, "Initial data has 0 likelihood.")
  }

  list(error = errortext, lik0 = lik0)
}