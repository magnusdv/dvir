# Data for Exercise 4.9.7 of Kling et al. (2021) 
# Mass Identifications: Statistical Methods in Forensic Genetics"

KETPex497 = familias2dvir("data-raw/KETPex497.fam", missingPrefix = "MP")
usethis::use_data(KETPex497, overwrite = TRUE)

# Check
if(FALSE){
  plotDVI(KETPex497)
  m = pairwiseLR(KETPex497)
  res = jointDVI(KETPex497)
  res
}