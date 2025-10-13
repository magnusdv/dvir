# Data for Exercise 4.9.7 of Kling et al. (2021) 
# Mass Identifications: Statistical Methods in Forensic Genetics"

KETPexer497 = familias2dvir("data-raw/KETPex497.fam", missingPrefix = "M")
usethis::use_data(KETPexer497, overwrite = TRUE)

# Check
if(FALSE){
  plotDVI(KETPexer497)
  dviSolve(KETPexer497)
}