# Dataset based on Exercise 3.3 EKM (2015)
# NB: OUTDATED CODE!
 
library(pedtools)

# Load fam file
x = readFam("data-raw/Solution3_3.fam", dedup = FALSE)

### PM data
pm = x$`Unidentified persons` |> 
  setSex(ids = founders, sex = 2) |>  # all female
  relabel(old = paste0("PM", 1:8), 
          new = paste0("V", 1:8))


### AM data
missing = paste0("M", 1:5)
refs = paste0("R", 1:5)

am = lapply(1:5, function(i) {
  fam = x[[paste("Family", i)]][["Reference pedigree"]] |> 
    relabel(old = "Missing person", new = missing[i])
  
  # Relabel typed member to R1, R2, ... (each family has just one)
  relabel(fam, old = typedMembers(fam), new = refs[i])
})

# Collect and save
planecrash = dviData(pm = pm, am = am, missing = missing)
usethis::use_data(planecrash, overwrite = TRUE)
