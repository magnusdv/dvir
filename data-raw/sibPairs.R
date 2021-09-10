library(forrel)
library(dvir)

pm = list()
victims = paste0("V", 1:20)
sex = rep(c(1,2), 10)

set1 = c("CSF1PO", "D2S1338", "D3S1358", "D5S818", "D7S820", "D8S1179", "D13S317", "D16S539",
         "D18S51", "D19S433", "D21S11", "FGA", "TH01", "TPOX", "VWA")
db = NorwegianFrequencies[set1]
db = NorwegianFrequencies
for (i in 1:20){
  pm[[i]] = singleton(victims[i],sex[i])
  pm[[i]]  = setMarkers(pm[[i]], locusAttributes = db)
}

am = list()
missing = paste0("M", 1:20)
x = addChildren(nuclearPed(1), father = 3, nch = 2, sex =c(1,2))

for (i in 1:10){
  am[[i]] = relabel(x, c((1:4)+4*(i-1), missing[2*i-1], missing[2*i]))
  am[[i]]  = setMarkers(am[[i]], locusAttributes = db)
}
gms = seq(2, 38, by = 4)
am = profileSim(am, ids = c(missing, gms), seed = 17)
pm = transferMarkers(am[[1]], pm, idsFrom = missing, idsTo = victims)
am = setAlleles(am[[1]], ids = missing, alleles = 0)

# Collect and save
sibPairs = list(pm = pm, am = am, missing = missing)

usethis::use_data(sibPairs, overwrite = TRUE)

# Examples
if(FALSE){
  # All but M10 correctly identified:
  res1 = sequentialDVI(pm, am, missing, updateLR = FALSE, verbose = TRUE)
  sequentialDVI(pm, am, missing, updateLR = TRUE, verbose = TRUE)
  # Takes a while, how long ? Reduce threshold?
  date();res3 = jointDVI(pm, am, missing, verbose = TRUE);date()
  n = ncomb(10, 10, 10, 10) 
  # 5.506636e+16 different assignments. Should be impossible without updates
}