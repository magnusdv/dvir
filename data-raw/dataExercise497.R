# Data for Exercise 4.9.7 of Kling et al. (2021) 
# Mass Identifications: Statistical Methods in Forensic Genetics"
# load data
rm(list = ls()) # removes all R objects
con <- url("http://familias.name/BookKETP/Files/dataExercise_4_9_7.RData") 
load(con) 
close(con) # Finished loading data: from, to, ids.to and moves
rm(con)

library(pedprobr)
library(forrel)
# Make data set with  equal mutation model with rate = 0.001 and rename variables
pm = list()
am = list()
missing = ids.to
r = 0.001
for (i in 1:3){
  pm[[i]] = setMutationModel(from[[i]], model = "equal", rate = r)
  am[[i]] = setMutationModel(to[[i]], model = "equal", rate = r)
}

# Check
 # plotPedList(list(am, pm), marker = 1, hatched = typedMembers, col = list(red = missing))

# Collect and save
dataExercise497 = list(pm = pm, am = am, missing = missing)

usethis::use_data(dataExercise497, overwrite = TRUE)
