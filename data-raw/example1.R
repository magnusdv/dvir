### DVI example showing that forward selection may fail

library(pedtools)

loc = list(name = "L1", 
           alleles = 1:3,
           afreq = c(0.05, 0.05, 0.9))

# PM data
victims = paste0("V", 1:3)

pm.df = data.frame(famid = victims, id = victims,
                   fid = 0, mid = 0, sex = 1,
                   L1 = c("1/1", "2/2", "1/2"))
pm = as.ped(pm.df, locusAttributes = loc)

# AM data
am = nuclearPed(father = "M1", mother = "R1", child = "M2")
am = addChildren(am, father = "M2", mother = "NN", ids = "M3")
L1 = marker(am, "R1" = "1/1", name = "L1", alleles = loc$alleles, afreq = loc$afreq)
am = setMarkers(am, L1)

missing = c("M1", "M2", "M3")

# Check
plotPedList(list(am, pm), marker = 1, hatched = typedMembers, col = list(red = missing))

# Collect and save
example1 = list(pm = pm, am = am, missing = missing)

usethis::use_data(example1, overwrite = TRUE)
