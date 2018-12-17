#' Find next candidates
#' 
#' New moves are generated, NULL
#' if none are possible, for one and two step
#' 
#' @param from A list of singeltons 
#' @param to A list of pedigrees 
#' @param ids.from Character vector with names of victims
#' @param ids.to Character vector with names of missing persons
#' @param one logic If TRUE, only one-step
#' @return list of matrices of one and two-moves
#' @examples 
#' ids.from = c("V1", "V2")
#' from = singletonList(2, ids = ids.from)
#' m = marker(from[[1]], alleles = 1:2, "V1" = 1)
#' from[[1]] = addMarkers(from[[1]], m)
#' m = marker(from[[2]], alleles = 1:2, "V2" = 1:2)
#' from[[2]] = addMarkers(from[[2]], m)
#' ids.to = c("MP1", "MP2")
#' to = nuclearPed(2, children = ids.to, father = "R1", mother = "R2")
#' m = marker(to, "R1" = 1:2, "R2"= 1:2)
#' to = addMarkers(to, m)
#' # The thre possibilities: a move of two, one or none:
#' tab = generate(from, to, ids.from, ids.to, one = FALSE)
#' from[[1]] = swapSex(from[[1]], "V1")
#' tab = generate(from, to, ids.from,ids.to)
#' from[[2]] = swapSex(from[[2]], "V2")
#' tab = generate(from, to, ids.from,ids.to)
#' @export
#' 
generate = function(from, to, ids.from, ids.to, one = FALSE){
  if(one)
    tab2  = NULL
  else
    tab2 = loopPair(from, to, ids.from, ids.to)
  amsex = data.frame(ids.to = ids.to, sex = getSex(to, ids.to))
  rownames(amsex) = ids.to
  pmsex = data.frame(ids.from = ids.from, sex = getSex(from, ids.from))
  rownames(pmsex) = ids.from
  sex = list(pmsex = pmsex, amsex = amsex)
  tab1 = loop2(ids.from,ids.to, sex)
  list(onestep = tab1, twostep = tab2)
}

#' @export
loopPair = function(from, to, ids.from, ids.to){
# Generate list of candidate pairs to move
# 
# All pairs of victims and pairs of missing persons
# that are possible accounting for sex are generated.
# ids.from Character vector with names of two victims
# ids.to Character vector with names of two missing persons
# from A list of singeltons 
# to A list of pedigrees 
  
  ids.from = sort(ids.from)
  ids.to = sort(ids.to)
  nf = length(ids.from)
  nt = length(ids.to)
   if (nf <= 1 | nt <= 1) return(NULL)
  n1 = arrangements::ncombinations(length(ids.from), 2)
  n2 = arrangements::ncombinations(length(ids.to), 2)
  comb1 = arrangements::combinations(ids.from, 2)
  foo1 = do.call("rbind", rep(list(comb1), n2))
  comb2 = arrangements::combinations(ids.to, 2)
  foo2 = matrix(rep(comb2, each = n1), ncol = 2)
  loop1 = cbind(foo1, foo2)
  # All unordered combinations have been generated
  # Find sex and discard impossible suggestions
  sex12 = apply(matrix(loop1[, 1:2], ncol = 2), 1, 
                function(x, from) getSex(from, x), from)
  sex12  = matrix(sex12, ncol = 2, byrow = TRUE)
  sex34 = apply(matrix(loop1[, 3:4], ncol = 2), 1, 
                function(x, to) getSex(to, x), to)
  sex34  = matrix(sex34, ncol = 2, byrow = TRUE)
  sex = cbind( sex12, sex34)
  index = apply(matrix(sex, ncol = 4), 1, 
                function(x) all(sort(x[1:2]) == sort(x[3:4])))
  if(sum(index) > 0)
    loop1 = matrix(cbind(loop1, sex)[index,], ncol = 8 )
  else
    return(NULL)
  # Add permutations when sex is the same and order output
  samesex = apply(matrix(loop1[,5:8], ncol = 4), 1, 
                  function(x) all(x == 1) | all(x ==2))
  if(sum(samesex) > 0){
    pe1 = c(1, 2, 4, 3, 5:8)
    pe2 = c(2, 1, 4, 3, 5:8)
    pe3 = c(2, 1, 3, 4, 5:8)
    combinations = rbind(loop1,loop1[samesex, pe1], loop1[samesex, pe2], 
                         loop1[samesex, pe3])
  } else
    combinations = loop1
  colnames(combinations) =  c("from1", "from2", "to1", "to2",
                           "sex1", "sex2", "sex3", "sex4")
  combinations = data.frame(combinations)
  combinations[do.call(order, combinations),]
  }

#' @export
loop2 <-
function(ids.from, ids.to, sex){
# Find single candidate list
    loop1 = loop(ids.from, ids.to)
    d1 = dim(loop1)[1]
    index = rep(NA, d1)
    for (i in 1:d1){
      from1 = as.character(loop1[i,1])
      to1 = as.character(loop1[i,2])
      index[i] = sex[[1]][from1, 2] == sex[[2]][to1, 2]
    }
    ni = sum(index)
    if (ni == 0) return(NULL)
    loop1 = matrix(loop1[index, ], ncol = 2)
    colnames(loop1) = c("from", "to")
    loop1
  }

loop <-
function(ids.from, ids.to){
    from2 = rep(ids.from, each = length(ids.to))
    to2 = rep(ids.to, length(ids.from))
    loop1 = cbind(from2, to2)
    loop1
  }
