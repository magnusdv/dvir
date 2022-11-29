# Data for Example 4.8.1 of Kling et al. (2021) 
# Mass Identifications: Statistical Methods in Forensic Genetics"

library(forrel)
x = readFam("data-raw/Ch4-table-4-4.fam", verbose = F)
V1 = x$H6[[1]]
V2 = x$H5[[2]]
pm = list(V1, V2)
am = x$H7[[1]]
am = relabel(am, c("MP2", "MP1", "R", "FA", "GM"))
am = setAlleles(am, ids = c("MP1", "MP2"), alleles = 0)
missing = paste0("MP", 1:2)

dataExample481 = dviData(pm, am, missing)

usethis::use_data(dataExample481, overwrite = TRUE)

# Check
if(FALSE) {
  plotDVI(dataExample481, marker = 1:2)
  res = jointDVI(dataExample481)
  res
}
