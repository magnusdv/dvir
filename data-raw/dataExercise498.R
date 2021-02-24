# Data for Exercise 4.9.8 of Kling et al. (2021) 
# Mass Identifications: Statistical Methods in Forensic Genetics"
# load data
rm(list = ls()) # removes all R objects
con <- url("http://familias.name/BookKETP/Files/dataExercise_4_9_8.RData") 
load(con) 
close(con) # Finished loading data: from, to, ids.to 
rm(con)
pm = from
am = to
missing = ids.to 

library(pedprobr)
library(forrel)
# Make data set with  equal mutation model with rate = 0.001 and rename variables
r = 0.001
for (i in 1:3)
  pm[[i]] = setMutationModel(pm[[i]], model = "equal", rate = r)
  am= setMutationModel(am, model = "equal", rate = r)


# Check
plotPedList(list(am, pm), marker = 1, hatched = typedMembers, col = list(red = missing))

# Collect and save
dataExercise498 = list(pm = pm, am = am, missing = missing)

usethis::use_data(dataExercise498, overwrite = TRUE)
