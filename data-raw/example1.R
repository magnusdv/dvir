### DVI example showing that forward selection may fail

library(pedtools)

loc = list(name = "L1", alleles = 1:3, afreq = c(0.05, 0.05, 0.9))

# PM data

pm = list(
  singleton("V1"),
  singleton("V2"),  
  singleton("V3")
) |> 
  addMarker(V1 = "1/1", V2 = "2/2", V3 = "1/2", locusAttr = loc)

# AM data

am = nuclearPed(father = "M1", mother = "R1", child = "M2") |> 
  addSon(parents = c("M2", "NN"), id = "M3") |> 
  addMarker(R1 = "1/1", locusAttr = loc)

missing = c("M1", "M2", "M3")

# Collect and save
example11 = dviData(pm = pm, am = am, missing = missing)
usethis::use_data(example1, overwrite = TRUE)

# Check
if(FALSE) {
  plotDVI(example1, marker = 1)
}
