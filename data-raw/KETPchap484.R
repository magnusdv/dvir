# Example 4.8.4 from Chapter 4 of Mass Identifications (KETP)

library(pedsuite)
x = pedFamilias::readFam("data-raw/Ch4-example-4-8-4.fam")

# Quick overview
plot(x, hatched = typedMembers, marker = 1)

vics    = paste0("V", 1:4)
missing = paste0("M", 1:4)

# Extract pm and am datasets
pm = extractSingletons(x, vics)

am = list(
  F1 = x[[1]], 
  F2 = mergePed(x[[2]], x[[4]], by = c("V2" = "added_1")), 
  F3 = x[[3]]) 

# Remove genotypes from victims and relabel them
am = am |> 
  removeGenotypes(ids = vics) |> 
  relabel(old = vics, new = missing) 

# Collect and save
KETPchap484 = dviData(pm = pm, am = am, missing = missing)

usethis::use_data(KETPchap484, overwrite = TRUE)

# Check
if(FALSE){
  plotDVI(KETPchap484, marker = 1:2)
  res = dviSolve(KETPchap484)
  res
}
