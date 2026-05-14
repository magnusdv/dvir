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

# More sensible ordering
ord = c("GF", "GM", "R1", "FA2", "MP1", "MP2", "R2", "MP3", "FA", 
        "R3", "GF2", "MO", "MP4", "MP5", "FA4", "R4", "UN2", "MO3", 
        "R5", "MP6", "MP7", "GM3", "MP8")
am = reorderPed(am, ord)

plot(am)

# Collect and save
grave = dviData(pm = pm, am = am, missing = missing)
usethis::use_data(grave, overwrite = TRUE)


# Check
if(FALSE){
  plotDVI(grave)
  dviSolve(grave)
}
