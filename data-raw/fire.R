# Example featured in the GLR paper (Egeland & Vigeland 2024)
# Code copied from https://github.com/magnusdv/GLR-paper/blob/master/examples/example2-dvi/1-generate-data.R

library(pedsuite)
library(dvir)

# Database
db = forrel::NorwegianFrequencies[1:23]

# Create empty DVI object
am = nuclearPed(ch = c("R1", "M1", 3), sex = c(1,2,1)) |>
  addChildren(3, ids = c("M2", "M3"), sex = 2) |>
  setMarkers(locusAttributes = db) |>
  setMutmod(model = "prop", rate = 0.001)

pm = singletons(c("V1", "V2", "V3"), sex = 2) |>
  transferMarkers(from = am, to = _)

dvi0 = dviData(pm, am, missing = c("M1", "M2", "M3"))

# True pairing
truth = c(V1 = "M1", V2 = "M2", V3 = "M3")

# Simulate profiles for R1 and victims
fire = dviSim(dvi0, refs = "R1", truth = truth, seed = 151)

usethis::use_data(fire, overwrite = TRUE)
