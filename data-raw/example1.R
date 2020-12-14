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
am = nuclearPed(father = "MP1", mother = "R1", child = "MP2")
am = addChildren(am, father = "MP2", mother = "NN", ids = "MP3", verbose = FALSE)
m = marker(am, "R1" = 1, alleles = als, afreq = p, name = "L1")
am = setMarkers(am, m)

missing = c("MP1", "MP2", "MP3")

# Check
plotPedList(list(am, pm), marker = 1, hatched = typedMembers, col = list(red = missing))

# Collect and save
example1 = list(pm = pm, am = am, missing = missing)

usethis::use_data(example1, overwrite = TRUE)
