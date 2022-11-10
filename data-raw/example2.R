### Small DVI example, 3 victims, two reference families with 3 MP-s

library(pedtools)

loc = list(name = "L1", 
           alleles = 1:10,
           afreq = rep(1/10,10))

# PM data
victims = paste0("V", 1:3)

pm.df = data.frame(famid = victims, id = victims,
                   fid = 0, mid = 0, sex = c(1, 1, 2),
                   L1 = c("1/1", "1/2", "3/4"))
pm = as.ped(pm.df, locusAttributes = loc)

# AM data
am1 = nuclearPed(father = "M1", mother = "R1", child = "M2")
L1 = marker(am1, "R1" = "2/2", name = "L1", alleles = loc$alleles, afreq = loc$afreq)
am1 = setMarkers(am1, L1)

am2 = nuclearPed(father = "R2", mother = "MO2", child = "M3", sex = 2)
L1 = marker(am2, "R2" = "3/3", name = "L1", alleles = loc$alleles, afreq = loc$afreq)
am2 = setMarkers(am2, L1)

am = list(am1, am2)

missing = c("M1", "M2", "M3")

# Check
plotPedList(list(am, pm), marker = 1, hatched = typedMembers, col = list(red = missing))

# Collect and save
example2 = dviData(pm = pm, am = am, missing = missing)

usethis::use_data(example2, overwrite = TRUE)
