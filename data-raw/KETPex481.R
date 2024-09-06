# Data for Example 4.8.1 of Kling et al. (2021) 
# Mass Identifications: Statistical Methods in Forensic Genetics"

library(dvir)
x = pedFamilias::readFam("data-raw/Ch4-table-4-4.fam", verbose = F)

# Missing
missing = paste0("MP", 1:2)

# PM data
pm = list(x$H6[[1]], x$H5[[2]])

# AM data
am = x$H7[[1]] |> 
  relabel(new = c("MP2", "MP1", "R1", "FA", "GM")) |> 
  setAlleles(ids = missing, alleles = 0)

KETPex481 = dviData(pm, am, missing)

usethis::use_data(KETPex481, overwrite = TRUE)

# Check
if(FALSE) {
  plotDVI(KETPex481, marker = 1:2)
  pairwiseLR(KETPex481)$LRmatrix
  dviSolve(KETPex481)
}
