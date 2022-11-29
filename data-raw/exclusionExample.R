# exclusionExample provided by Lourdes Prieto

exclusionExample = dvir::familias2dvir("exclusionExample.fam",
                                  victimPrefix ="V",
                                  missingPrefix = "M",
                                  refPrefix = "R")

library(pedtools)
usethis::use_data(exclusionExample, overwrite = TRUE)


# Check
if(FALSE){
  plotDVI(exclusionExample, am = 1)
  m = pairwiseLR(exclusionExample)
  res = jointDVI(exclusionExample)
  res
}