#' Generate list of candidate pairs to move
#' 
#' All pairs of victims and pairs of missing persons
#' that are possible accounting for sex are generated.
#' @param from1 Character vector with names of victims
#' @param to1 Character vector with names of missing persons
#' @param pm A list of singeltons 
#' @param am A list of pedigrees 
#' @return A data frame with eight columns. The first two describe a candidate pair
#' to bemoved. The next two a destination pair. Then follows sex information. 
#' @details The output is ordered starting from the first column.
#' NULL is returned if there are no possible moves.
#' @export
#' @examples
#' from1 = paste("V", 1:3, sep = "")
#' pm = singletonList(3, sex = c(1, 1, 2))
#' to1 = paste("MP", 1:4, sep = "")
#' am = nuclearPed(4, children = to1, sex = c(1, 1, 2, 2))
#' res = loopPair(from1, to1, pm, am)
#' 
loopPair = function(from1, to1, pm, am){
  from1 = sort(from1)
  to = sort(to1)
  nf = length(from1)
  nt = length(to1)
   if (nf <= 1 | nt <= 1) return(NULL)
  n1 = arrangements::ncombinations(length(from1), 2)
  n2 = arrangements::ncombinations(length(to1), 2)
  comb1 = arrangements::combinations(from1, 2)
  foo1 = do.call("rbind", rep(list(comb1), n2))
  comb2 = arrangements::combinations(to1, 2)
  foo2 = matrix(rep(comb2, each = n1), ncol = 2)
  loop1 = cbind(foo1, foo2)
  # All unordered combinations have been generated
  # Find sex and discard impossible suggestions
  sex12 = apply(matrix(loop1[, 1:2], ncol = 2), 1, 
                function(x, pm) getSex(pm,x), pm)
  sex12  = matrix(sex12, ncol = 2, byrow = TRUE)
  sex34 = apply(matrix(loop1[, 3:4], ncol = 2), 1, 
                function(x, am) getSex(am,x), am)
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
