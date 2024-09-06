# Family grave data in Kling et al. (2021) 
# "Mass Identifications: Statistical Methods in Forensic Genetics"

library(dvir)
load(url("http://familias.name/BookKETP/Files/Grave.RData"))

pm = from
am = to
missing = ids.to

# Relabel
old = labels(am)[c(1:5,7:8,16,20:21)]
new = c("MO", "FA", "GM", "GF", "GF2", "MO3", "FA2", "UN2", "GM3", "FA4")
am = relabel(am, new, old)
plot(am)

# Collect and save
grave = dviData(pm = pm, am = am, missing = missing)
usethis::use_data(grave, overwrite = TRUE)


# Check
if(FALSE){
  plotDVI(grave)
  dviSolve(grave)
}