x = forrel::readFam("http://familias.name/BookKETP/Files/Ch4-example-4-8-4.fam")
V1 = branch(x$H1[[1]],"V1")
V2 = branch(x$H1[[2]],"V2")
V3 = branch(x$H1[[4]],"V3")
V4 = branch(x$H1[[3]],"V4")
pm = list(V1, V2, V3, V4)

R1 = branch(x$H1[[1]],"R1")
R2 = branch(x$H1[[2]],"R2")
R3 = branch(x$H1[[3]],"R3")

am1 = nuclearPed(father = "R1", mother  = "MO1", child = "MP1", sex = 2)
am1 = transferMarkers(R1, am1)
am2 = nuclearPed(father = "R2", mother  = "GM2", child = "MP2")
am2 = addChildren(am2, father ="MP2", mother = "MO2", ids= "MP3")
am2 = transferMarkers(R2, am2)

am3 = nuclearPed(father = "GF3", mother  = "GM3", children = c("FA3", "MP4"), sex = c(1,2))
am3 = addChildren(am3, father ="FA3", mother = "MO3", ids= "R3")
am3 = transferMarkers(R3, am3)
am = list(am1, am2, am3)
missing = paste0("MP", 1:4)

# Collect and save
dataCh4 = dviData(pm = pm, am = am, missing = missing)

usethis::use_data(dataCh4, overwrite = TRUE)

# Check
if(FALSE){
  plotDVI(dataCh4, marker = 1:2)
  res = jointDVI(dataCh4, disableMutations = FALSE)
  res[c(1,2,30,49),]
}
