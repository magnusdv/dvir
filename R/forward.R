#' Forward DVI selection
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
#' @param from A list of singeltons 
#' @param to A list of pedigrees
#' @param ids.from Character vector with names of victims
#' @param ids.to Character vector with names of missing persons
#' @param LRlimit Double. Threshold for LR
#' @param eliminate logic If TRUE, irrelevant persons are removed
#' @param singleStep logic If TRUE, one step analysis
#' @param twoStep logic If TRUE, two step analysis
#' @return List with three elements corresponding to the three items
#' described above 
#' @export
#' @examples
#' 
#' als = 1:2
#' freq = c(0.5, 0.5)
#' mutmodel = "proportional"
#' rate = 0.00
#' v1 = singleton("V1")
#' m = marker(v1, alleles = als, afreq = freq, "V1" = c(1, 1),
#' mutmod = mutmodel, rate = rate)
#' v1 = addMarkers(v1, m)
#' v2 = singleton("V2")
#' m = marker(v2, alleles = als, afreq = freq, "V2" = c(1, 2),
#' mutmod = mutmodel, rate = rate)
#' v2 = addMarkers(v2, m)
#' v3 = singleton("V3", 2)
#' m = marker(v3, alleles = als, afreq = freq, "V3" = c(1, 2),
#'           mutmod = mutmodel, rate = rate)
#' v3 = addMarkers(v3, m)
#' from = list(v1, v2, v3)
#' to = nuclearPed(3, sex = c(1,1,2), father = "fa", mother = "mo", children = c("MP1", "MP2", "MP3"))
#' m = marker(to, alleles = als, afreq = freq, "fa" = 1:2, "mo" = c(1,1),
#'         mutmod = mutmodel, rate = rate)
#' to = addMarkers(to, m)
#' plotPedList(list(from, to), marker = 1)
#' ids.from = c("V1", "V2", "V3")
#' ids.to = c("MP1","MP2", "MP3")
#' res = forward(from, to, ids.from, ids.to, LRlimit = -1)

#' from = singletonList(4)
#' ids.from = unlist(lapply(from, function(x) x$ID))
#' ids.to = c("MP1", "MP2", "MP3")
#' to = nuclearPed(3, children = ids.to, father = "R1", mother = "R2")
#' m = marker(to, "R1" = 1, "R2"= 1, alleles = 1:2)
#' to = addMarkers(to, m)
#' limit = -1
#' #' plotPedList(list(from, to), marker = 1)
#' res = forward(from, to, ids.from, ids.to, LRlimit = -1, 
#'      eliminate = TRUE,singleStep = FALSE, twoStep = FALSE) 
#' res = forward(from, to, ids.from, ids.to, LRlimit = -1, 
#'      eliminate = TRUE,singleStep = FALSE, twoStep = TRUE) 
#' res = forward(from, to, ids.from, ids.to, LRlimit = -1, 
#'      eliminate = FALSE, singleStep = FALSE, twoStep = TRUE) 

forward <-
function(from, to, ids.from , ids.to, LRlimit = 1, eliminate = TRUE, 
         singleStep = TRUE, twoStep = TRUE){
  
  index1 = function(x, from){
    x[1] == from[1] | x[1] == from[2] | x[2] == from[1] | x[2] == from[2]
  }
  # Check input and initialise for output
  check = checkInput(from, to, ids.from, ids.to)
  if(!is.null(check$error)) stop(check$error)
  lik0 = check$lik0
  finished = FALSE
  eliminateTable = NULL # Table from reduce, i.e. elimation
  oneOnly = NULL # One step analysis output
  one = NULL #One step after two step
  summaryTable = NULL # Final summary table if twoStep

  #Unconditional matching. Remove irrelevant victims and missing persons
  if (eliminate){ 
    eliminateTable = reduce(from, to, ids.from, ids.to, 
                          lik0 = lik0, limit = LRlimit)
    if(is.null(eliminateTable))
      stop("No identifications possible")
    #Finished if only one identification possible
    ids.from = as.character(unique(eliminateTable[,1]))
    ids.to = as.character(unique(eliminateTable[,2]))
    finished = dim(eliminateTable)[1] < 2
    }
  # One-at-time identification. finish = TRUE if two or less identified
  if(singleStep & !finished){
    oneOnly = oneStep(from, to, ids.from , ids.to, 
                      LRlimit = LRlimit, lik0 = lik0)
    finished = dim(oneOnly[[1]])[1] <= 3
  }
  if(!twoStep | finished){
    ret = list(eliminateTable = eliminateTable,
               oneStepResults = oneOnly, summaryTable = summaryTable, 
               oneFinish = one)
    return(ret)
  }
  #Initialise
  tab2 = loopPair(from, to,ids.from, ids.to)
  finished = ifelse(is.null(tab2), TRUE, FALSE)
  if(!is.null(tab2)) tab2 = as.matrix(tab2[,1:4], ncol = 4)
  mi = min(length(ids.from), length(ids.to))
  i = 1
  fromList = NA
  toList = NA
  likList = check$lik0
  while(!finished & i <= mi){
      two = apply(matrix(tab2, ncol = 4), 1, 
            function(x, from, to) lik1(from, to, x[1:2], x[3:4]), from, to)
      likresTwo = unlist(lapply(two, function(x) x$lik))
      noTwo = as.integer(which.max(likresTwo))
      likTwo = likresTwo[[noTwo]]
      #update 
      if(likTwo == 0)
        finished = TRUE
      else{
        from1 = tab2[noTwo, 1:2]
        fromList = c(fromList, from1)
        to1 = tab2[noTwo, 3:4]
        toList = c(toList, to1)
        likFirst = lik1(from, to, from1[1], to1[1])$lik
        likList = c(likList, c(likFirst, likTwo))
        ind1 = apply(matrix(tab2[, 1:2], ncol = 2), 1, index1, from1)
        ind2 = apply(matrix(tab2[, 3:4], ncol=2), 1, index1, to1)
        ind3 = !ind2 & !ind1
        ids.from = setdiff(ids.from, from1)
        ids.to = setdiff(ids.to, to1)
        from = two[[noTwo]]$pm2
        to = two[[noTwo]]$am2
        if(sum(ind3) > 0){
            tab2 = matrix(tab2[ind3, ], ncol = 4)
            i = i+1
        } else
            finished =TRUE
        }
  }
   d1 = 1
   check = checkInput(from, to, ids.from, ids.to)
   if(is.null(check$error)){
      one = oneStep(from, to, ids.from , ids.to, LRlimit = LRlimit, 
                    lik0 = check$lik0, maxid = NULL)
      oneTab = one[[1]]
      d1 = dim(oneTab)[1]
   } else
     one = NULL

   if(!is.null(one) & d1 > 1){
     oneTab = one[[1]]
     fromCol = c(fromList, levels(oneTab[,1]))
     toCol = c(toList, levels(oneTab[,2]))
     d1 = dim(oneTab)[1]
     if (d1 > 1)
        likCol = c(likList, oneTab[-1,3])
     else
       likCol = likList
     steps = c(0, rep(1:2, length(likList[-1])/2))
     steps = c(steps, rep(1, d1-1))
     summaryTable = data.frame(from = fromCol, to = toCol, lik = likCol,
                               steps = steps)
   } else {
     summaryTable = data.frame(from = fromList, to = toList, lik = likList,
                               steps = c(0,rep(1:2, length(likList[-1])/2) ))
   }
   rownames(summaryTable) = 0:(dim(summaryTable)[1]-1)  
   ret = list(eliminateTable = eliminateTable,
              oneStepResults = oneOnly, summaryTable = summaryTable, 
              oneFinish = one)
   ret
}


 