### DVI case illustrating symmetric solutions

library(pedsuite)

refs = paste0("R", 1:3)
miss = paste0("M", 1:6)
vics = paste0("V", 1:6)

fams = list(
  F1 = nuclearPed(fa = "R1", mo = 1, children = c("M1", "M2")),
  F2 = nuclearPed(fa = "R2", mo = 2, children = c("M3", "M4")),
  F3 = nuclearPed(fa = "R3", mo = 3, children = c("M5", "M6"))
)

# Database
db = NorwegianFrequencies[1:13]

# Simulate AM data
am = fams |> profileSim(ids = c(refs, miss), markers = db, seed = 1)

# Transfer to PM data
pm = singletons(miss) |> 
  transferMarkers(from = am, to = _) |> 
  relabel(new = vics)

# Remove data from missing
am = am |> setAlleles(ids = miss, alleles = 0)


### Modify to illustrate different kinds of undecidability --------------------

# Remove victim V4 (Idea: Now V3 matches both M3 and M4)
pm2 = pm |> removeIndividuals("V4")

# Remove M6 (Idea: Now V5 and V6 both match M5)
am2 = am |> removeIndividuals("M6")
miss2 = setdiff(miss, "M6")

# DVI object
symmetricSibs = dviData(pm = pm2, am = am2, missing = miss2)
usethis::use_data(symmetricSibs, overwrite = TRUE)

if(F) {
  dvi = symmetricSibs
  plotDVI(dvi)
  u = findUndisputed(dvi)
  e = findExcluded(dvi)
  dviJoint(e$dviReduced)
  amDrivenDVI(dvi)
}
