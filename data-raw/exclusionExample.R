# exclusionExample provided by Lourdes Prieto

rm(list = ls()) # removes all R objects
con <- url("https://familias.name/dviapp/exclusionExample.RData") 
load(con) 
close(con) # Finished loading data: from, to, ids.to and moves
rm(con)

library(pedtools)

# Collect and save
exclusionExample = dviData(pm = exclusionExample$pm, 
                           am = exclusionExample$am, 
                           missing = exclusionExample$missing)
usethis::use_data(exclusionExample, overwrite = TRUE)


# Check
if(FALSE){
  plotDVI(exclusionExample, am = 1)
  m = pairwiseLR(exclusionExample)
  res = jointDVI(exclusionExample)
  res
}