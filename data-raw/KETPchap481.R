# Data for Example 4.8.1 of Kling et al. (2021) 
# Mass Identifications: Statistical Methods in Forensic Genetics"

x = pedFamilias::readFam("data-raw/Ch4-table-4-4.fam")

plotPedList(x, dev.height = 5, margins = 1, hatched = typedMembers, marker = 1:2)

vics    = paste0("V", 1:2)
missing = paste0("M", 1:2)

# PM data
pm = extractSingletons(x$H1, vics)

# AM data
am = x$H7[[1]] |> 
  removeGenotypes(ids = vics) |> 
  parentsBeforeChildren() |> 
  relabel(c(V1 = "M2", 
            V2 = "M1",
            R = "R1",
            added_1 = "FA", 
            added_2 = "GM")) 
  

# Collect and save
KETPchap481 = dviData(pm, am, missing)

usethis::use_data(KETPchap481, overwrite = TRUE)

# Check
if(FALSE) {
  plotDVI(KETPchap481, marker = 1:2)
  pairwiseLR(KETPchap481)$LRmatrix
  dviSolve(KETPchap481)
  jointDVI(KETPchap481, verb = F)
}
