#' Exercise 3.3 of EKM
#' 
#' The exercises 3.3 of EKM is studied vie simulation
#' 
#' @param from A list of singeltons 
#' @param to A list of pedigrees
#' @param ids.from Character vector with names of victims
#' @param ids.to Character vector with names of missing persons
#' @param simID A list of mp-s to simulate
#' @param moveID A list of victims to be replaced by the simulkated persons
#' @param N integer No of simulations
#' @param limit Double
#' @param seed Integer
#' @param mutRateSim double
#' @param mutRateCalc double
#' 
#' @details The simulated results are evaluated by the simple marginal search which
#' is sufficient for this example because of the informative realtives
#' @export
#' @examples
#' 
#' data(dvi.3.3.nomut) # Without mutations
#' dat = dvi.3.3.nomut
#' from = dat$pm
#' to = dat$am
#' ids.from = dat$vict
#' ids.to =  dat$miss
#' simID = ids.to
#' moveID = ids.from[1:5]
#' x1 = example.3.3(from, to, ids.from, ids.to, simID, moveID, 
#' N = 100, seed = 123, limit = 1,mutRateSim = 0.003, mutRateCalc = 0)
#' 
#' x2 = example.3.3(from, to, ids.from, ids.to, simID, moveID, 
#' N = 100, seed = 123, limit = 1,mutRateSim = 0.003, mutRateCalc = 0.003)
#' 
#' Plot
#' windows()
#' plotPedList(c(to), newdev = FALSE, marker=1, skip.empty.genotypes = TRUE, frames = T, 
#' frametitles = paste("F",1:5, sep="") )
#' plotPedList(c(from), newdev = FALSE, marker=1, 
#' skip.empty.genotypes = TRUE,  frames = 1:8)

example.3.3 = function(from, to, ids.from, ids.to, simID, moveID, N = 1, 
                       limit = 1, seed = NULL, mutRateSim = 0, mutRateCalc = 0){
  f = function(x){
    x1 = split(x, x$from)
    x2 = lapply(x1, function(z) z[1,])
    x3 = lapply(x2, function(z) {
      a = z[1]
      b = z[2]
      if(a == "PM1" & b == "MP1")   z[1] = "MP1"
      else if(a == "PM2" & b == "MP2") z[1] = "MP2"
      else if(a == "PM3" & b == "MP3") z[1] = "MP3"
      else if(a == "PM4" & b == "MP4") z[1] = "MP4"
      else if(a == "PM5" & b == "MP5") z[1] = "MP5"
      else z[1:2] = NA
      z
    })
    sum(unlist(lapply(x3, function(z) z[1]==z[2])), na.rm =T)== 5
  }
  
  if(length(simID) != length(moveID))
    stop("length simID and moveID must be equal")
  for (i in 1:5)
    for (j in 1:15) 
      mutmod(to[[i]]$markerdata[[j]]) =  mutationModel("proportional", 
                                alleles = alleles(to[[i]]$markerdata[[j]]),
                                afreq = afreq(to[[i]]$markerdata[[j]]), rate = mutRateSim)
  from = setAlleles(from, alleles = 0, ids = ids.from)
  sim = profileSim(c(from, to), N = N, 
                   ids = c(ids.from, simID), conditions = 1:15, seed = seed)
  tab1 = list()
  for (i in 1:N){
    g = getAlleles(sim[[i]], ids = simID)
    rownames(g) = moveID
    sim[[i]] = setAlleles(sim[[i]],ids =  moveID, alleles = g)
    sim[[i]] = setAlleles(sim[[i]],ids =  simID, alleles = 0)
    for (l in 1:13)
      for (s in 1:15) 
        mutmod(sim[[i]][[l]]$markerdata[[s]]) =  mutationModel("proportional", 
                                alleles = alleles(sim[[i]][[l]]$markerdata[[s]]),
                                afreq = afreq(sim[[i]][[l]]$markerdata[[s]]), rate = mutRateCalc)
    tab1[[i]] = reduce(sim[[i]][1:8], sim[[i]][9:13], ids.from, ids.to, limit = limit)
  }
  list(fracCorrect = sum(unlist(lapply(tab1, f))/N), tabs = tab1)
}

