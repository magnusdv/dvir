# Data for Exercise 4.9.8 of Kling et al. (2021) 
# Mass Identifications: Statistical Methods in Forensic Genetics"
# load data
rm(list = ls()) # removes all R objects
con <- url("http://familias.name/BookKETP/Files/dataExercise_4_9_8.RData") 
load(con) 
close(con) # Finished loading data: from, to, ids.to 
rm(con)

# Make data set with  equal mutation model with rate = 0.001 and rename variables
library(pedprobr)
r = 0.001
pm = setMutationModel(from, model = "equal", rate = r)
am = setMutationModel(to, model = "equal", rate = r)
missing = ids.to 

# Collect and save
dataExercise498 = dviData(pm = pm, am = am, missing = missing)
usethis::use_data(dataExercise498, overwrite = TRUE)


# Check
# plotPedList(list(am, pm), marker = 1, hatched = typedMembers, col = list(red = missing))
