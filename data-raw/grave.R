# Family grave data in Kling et al. (2021) 
# "Mass Identifications: Statistical Methods in Forensic Genetics"

rm(list = ls()) # removes all R objects
con <- url("http://familias.name/BookKETP/Files/Grave.RData") 
load(con) 
close(con) # Finished loading data: from, to, ids.to and moves
rm(con)

library(pedtools)
pm = from
am = to
missing = ids.to
old = labels(am)
new = old
new[c(1:5,7:8,16,20:21)] = c("MO", "FA", "GM", "GF", "GF2", "MO3", "FA2", "UN2", "GM3", "FA4")
am = relabel(am, new, old)

# Collect and save
grave = dviData(pm = pm, am = am, missing = missing)
usethis::use_data(grave, overwrite = TRUE)


# Check
if(FALSE){
  plotDVI(grave)
  m = pairwiseLR(grave)
  res = jointDVI(grave)
  res
}