# Dataset based on Exercise 3.3 EKM (2015)

library(pedtools)

# Load fam file
x = forrel::readFam("data-raw/Solution3_3.fam")

### PM data
pm = x$`Unidentified persons`

# Rename singletons: victims
victims = paste0("V", 1:8)
names(pm) = victims
pm = relabel(pm, old =  paste0("PM", 1:8), new = victims)

# Change sex of first 6 victims, so that all 8 are female
pm[1:6] = lapply(pm[1:6], function(v) swapSex(v, labels(v)))

### AM data
missing = paste0("M", 1:5)
refs = paste0("R", 1:5)

am = lapply(1:5, function(i) {
  fam = x[[paste("Family", i)]][["Reference pedigree"]]
  
  # Relabel typed member to R1, R2, ... (each family has just one)
  fam = relabel(fam, old = typedMembers(fam), new = refs[i])
  
  # Relabel missing person
  relabel(fam, old = "Missing person", new = missing[i])
})

# Collect and save
planecrash = list(pm = pm, am = am, missing = missing)

usethis::use_data(planecrash, overwrite = TRUE)
