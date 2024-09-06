# Dataset based on Exercise 3.3 EKM (2015)
# NB: OUTDATED CODE!
 
library(dvir)

# Load fam file
x = pedFamilias::readFam("data-raw/Solution3_3.fam", prefixAdded = "A")

missing = paste0("M", 1:5)
refs = paste0("R", 1:5)
vics = paste0("V", 1:8)

### PM data
pm = x$`Unidentified persons` |> 
  setSex(ids = founders, sex = 2) |>  # all female
  relabel(new = vics)


### AM data
am = lapply(1:5, function(i) {
  fam = x[[paste("Family", i)]] |> 
    relabel(old = "Missing person", new = missing[i]) 
  # refs -> R1, R2, ... (single ref in each fam)
  relabel(fam, old = typedMembers(fam), new = refs[i])
})

# Collect and save
planecrash = dviData(pm = pm, am = am, missing = missing)
checkDVI(planecrash)
usethis::use_data(planecrash, overwrite = TRUE)

if(FALSE) {
  plotDVI(planecrash)
  dviSolve(planecrash)
}