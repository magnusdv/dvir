# Data for Example 4.8.3 of Kling et al. (2021) 
# Mass Identifications: Statistical Methods in Forensic Genetics"
# load data
library(pedtools)
library(pedprobr)
library(forrel)
load("C:\\Users\\theg\\Dropbox\\MyBookKEP\\ch4\\old\\icmp2.Rdata")
pm = icmp2$from
am = icmp2$to
missing = paste("MP", 1:12, sep ="")


# Make data set with  equal mutation model with rate = 0.001 and rename variables
r = 0.001
for (i in 1:5)
  pm[[i]] = setMutationModel(pm[[i]], model = "equal", rate = r)
am = setMutationModel(am, model = "equal", rate = r)


# Check
library(dvir)
plotPedList(list(am, pm), marker = 1:2, hatched = typedMembers, col = list(red = missing))
plotPedList(list(am), marker = 1:2, hatched = typedMembers, col = list(red = missing))
# res = jointDVI(pm, am, missing)
# write.table(head(res[,-6]), file = "Table-4-5.txt", quote = F)
# Collect and save
dataExample483 = list(pm = pm, am = am, missing = missing)

usethis::use_data(dataExample483, overwrite = TRUE)
