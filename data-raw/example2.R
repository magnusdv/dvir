### Small DVI example, 3 victims, two reference families with 3 MP-s

library(pedtools)

# Attributes of the marker locus
loc = list(name = "L1", alleles = 1:10, afreq = rep(0.1, 10))

# PM data

pm = list(
  singleton("V1", sex = 1),
  singleton("V2", sex = 1),  
  singleton("V3", sex = 2)
) |> 
  addMarker(V1 = "1/1", V2 = "1/2", V3 = "3/4", locusAttr = loc)

# AM data

am = list(
  nuclearPed(father = "M1", mother = "R1",  child = "M2"),
  nuclearPed(father = "R2", mother = "MO2", child = "M3", sex = 2)
) |> 
  addMarker(R1 = "2/2", R2 = "3/3", locusAttr = loc)

# Missing persons
missing = c("M1", "M2", "M3")

# Create and save DVI object
example2 = dviData(pm = pm, am = am, missing = missing)
usethis::use_data(example2, overwrite = TRUE)

# Check
if(FALSE) {
  plotDVI(example2, marker = 1)
}
