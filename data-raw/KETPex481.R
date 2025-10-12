# Data for Example 4.8.1 of Kling et al. (2021) 
# Mass Identifications: Statistical Methods in Forensic Genetics"

x = pedFamilias::readFam("data-raw/Ch4-table-4-4.fam")

plotPedList(x, dev.height = 5, margins = 1, hatched = typedMembers, marker = 1:2)

vics    = paste0("V",  1:2)
missing = paste0("MP", 1:2)

# PM data
pm = extractSingletons(x$H1, vics)

# AM data
am = x$H7[[1]] |> 
  removeGenotypes(ids = vics) |> 
  parentsBeforeChildren() |> 
  relabel(c(V1 = "MP2", 
            V2 = "MP1",
            R = "R1",
            added_1 = "FA", 
            added_2 = "GM")) 
  

# Collect and save
KETPex481 = dviData(pm, am, missing)

usethis::use_data(KETPex481, overwrite = TRUE)

# Check
if(FALSE) {
  plotDVI(KETPex481, marker = 1:2)
  pairwiseLR(KETPex481)$LRmatrix
  dviSolve(KETPex481)
}
