# Data for Exercise 4.9.8 of Kling et al. (2021) Mass Identifications:
# Statistical Methods in Forensic Genetics" Read data and change sex, names to
# agree with previous dataExercise498

KETPex498 = familias2dvir("KETPex498.fam", missingPrefix ="MP")
KETPex498$am[[1]] = relabel(KETPex498$am[[1]], c("MP1", "MP2", "MP3"),
                            c("MP3", "MP1", "MP2"))
KETPex498$am[[1]] = swapSex(KETPex498$am[[1]], c("MP2", "MP3"))


# Collect and save
KETPex498 = dviData(pm = KETPex498$pm, 
                    am = KETPex498$am, 
                    missing = KETPex498$missing)
usethis::use_data(KETPex498, overwrite = TRUE)

# Check
if(FALSE){
  plotDVI(KETPex498, nrowPM = 3)
  m = pairwiseLR(KETPex498)
  res = jointDVI(KETPex498)
  res
}


