#' Forward DVI selection, main function
#' 
#' The input is checked and the NULL likelihood calculated.
#' Three steps can be included depending on the choice of 
#' the corresponding the three logic variables
#' 1) Impossible victims and missing persons can be removed based on 
#' genotypes (if mutations are not modelled) or sex.
#' 2) A one step selection can be performed in which case an identification
#' is made if the LR threshold is met. The data is updated.
#' 3) Pairwise identifications are done as long as possible and followed by
#' one steps if required.
#' 
#' @param pm A list of singeltons 
#' @param am A list of pedigrees
#' @param vp Character vector with names of victims
#' @param mp Character vector with names of missing persons
#' @param LRlimit Double. Threshold for LR
#' @param reduser logic If TRUE, irrelevant persons are removed
#' @param singleStep logic If TRUE, one step analysis
#' @param twoStep logic If TRUE, two step analysis
#' @return List with three elements corresponing to the three items
#' described above 
#' @examples
#' pm = singletonList(4)
#' vp = unlist(lapply(pm, function(x) x$ID))
#' m = marker(pm[[1]], alleles = 1:2, "V1" = 1)
#' pm[[1]] = addMarkers(pm[[1]], m)
#' m = marker(pm[[2]], alleles = 1:2, "V2" = 1)
#' pm[[2]] = addMarkers(pm[[2]], m)
#' m = marker(pm[[3]], alleles = 1:2, "V3" = 2)
#' pm[[3]] = addMarkers(pm[[3]], m)
#' m = marker(pm[[4]], alleles = 1:2, "V4" = 1)
#' pm[[4]] = addMarkers(pm[[4]], m)
#' mp = c("MP1", "MP2", "MP3")
#' am = nuclearPed(3, children = mp, father = "R1", mother = "R2")
#' m = marker(am, "R1" = 1, "R2"= 1, alleles = 1:2)
#' am = addMarkers(am, m)
#' limit = -1
#' #' plotPedList(list(pm, am), marker = 1)
#' res = forward(pm, am, vp, mp, LRlimit = -1, 
#'      reduser = TRUE,singleStep = FALSE, twoStep = FALSE) 
#' res = forward(pm, am, vp, mp, LRlimit = -1, 
#'      reduser = TRUE,singleStep = FALSE, twoStep = TRUE) 
#' res = forward(pm, am, vp, mp, LRlimit = -1, 
#'      reduser = FALSE, singleStep = FALSE, twoStep = TRUE) 
#' @export

forward <-
function(pm, am, vp , mp, LRlimit = 1, reduser = FALSE, singleStep = FALSE, twoStep = TRUE){
  # Check input and initialise for output
  check = checkInput(pm, am, vp, mp)
  lik0 = check$lik0
  if(!is.null(check$error)) stop(check$error)
  finished = FALSE
  unconditionalTable = NULL
  summaryTable = NULL
  one = NULL
  if (reduser){ #Remove irrelevant victims and missing persons
    res = reduce(pm, am, vp, mp, likNULL = lik0)
    if(is.null(res))
      stop("No identifications possible")
    #Finished if only one identification possible
    vp = res$vp
    mp = res$mp
    unconditionalTable = res
    finished = dim(unconditionalTable$ret)[1] < 2
    }
  # One-at-time identification and finish if no identified one or two
  if(singleStep & !finished){
    one = oneStep(pm, am, vp , mp, LRlimit = 0, lik0 = lik0)
    finished = dim(one[[1]])[1] <= 3
  }
  if(!twoStep | finished){
    ret = list(unconditionalTable = unconditionalTable,
               oneStepResults = one, summaryTable = summaryTable)
    return(ret)
  }
  #Initialise
  tab2 = loopPair(vp, mp, pm, am)
  finished = ifelse(is.null(tab2), TRUE, FALSE)
  if(!is.null(tab2)) tab2 = as.matrix(tab2[,1:4], ncol = 4)
  mi = min(length(vp), length(mp))
  i = 1
  fromList = NA
  toList = NA
  likList = check$lik0
  while(!finished & i <= mi){
      two = apply(tab2, 1, 
            function(x, pm, am) lik1(x[1:2], 
            x[3:4], pm, am), pm, am)
      likresTwo = unlist(lapply(two, function(x) x$lik))
      noTwo = as.integer(which.max(likresTwo))
      likTwo = likresTwo[[noTwo]]
      #update 
      if(likTwo == 0)
        finished = TRUE
      else{
        from = tab2[noTwo, 1:2]
        fromList = c(fromList, from)
        to = tab2[noTwo, 3:4]
        toList = c(toList, to)
        likFirst = lik1(from[1], to[1], pm, am)$lik
        likList = c(likList, c(likFirst, likTwo))
        ind1 = apply(tab2[,1:2], 1, index1, from)
        ind2 = apply(tab2[,3:4], 1, index1, to)
        ind3 = !ind2 & !ind1
        vp = setdiff(vp, from)
        mp = setdiff(mp, to)
        pm = two[[noTwo]]$pm2
        am = two[[noTwo]]$am2
        if(sum(ind3) > 0){
            tab2 = tab2[ind3, ]
            i = i+1
        } else
            finished =TRUE
        }
  }
   check = checkInput(pm, am, vp, mp)
   if(is.null(check$error)){
      one = oneStep(pm, am, vp , mp, LRlimit = 1, 
                    lik0 = check$lik0, maxid = NULL)
      oneTab = one[[1]]
      d1 = dim(oneTab)[1]
   } else
     one = NULL
   
   if(!is.null(one)){
     oneTab = one[[1]]
     fromCol = c(fromList, levels(oneTab[,1]))
     toCol = c(toList, levels(oneTab[,2]))
     d1 = dim(oneTab)[1]
     if (d1 > 1)
        likCol = c(likList, oneTab[-1,3])
     else
       likCol = likList
     steps = c(0, rep(1:2, length(likList[-1])/2))
     if(!is.null(d1) & d1>1) steps = c(steps, rep(1, d1-1))
     summaryTable = data.frame(from = fromCol, to = toCol, lik = likCol,
                               steps = steps)
   } else {
     summaryTable = data.frame(from = fromList, to = toList, lik = likList,
                               steps = c(0,rep(1:2, length(likList[-1])/2) ))
   }
   rownames(summaryTable) = 0:(dim(summaryTable)[1]-1)  
   ret = list(unconditionalTable = unconditionalTable, 
              oneStepResults = one, summaryTable = summaryTable)
   ret
}

 index1 = function(x, from){
   x[1] == from[1] | x[1] == from[2] | x[2] == from[1] | x[2] == from[2]
 }
 