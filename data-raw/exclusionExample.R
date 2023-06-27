# exclusionExample provided by Lourdes Prieto

exclusionExample = dvir::familias2dvir("data-raw/exclusionExample.fam",
                                        victimPrefix ="V",
                                        missingPrefix = "M",
                                        refPrefix = "R")

usethis::use_data(exclusionExample, overwrite = TRUE)


# Check
if(FALSE){
  plotDVI(exclusionExample, am = 1)
  m = pairwiseLR(exclusionExample)
  e = findExcluded(exclusionExample)
  res = jointDVI(exclusionExample)
  res
}