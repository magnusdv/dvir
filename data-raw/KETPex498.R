# Data for Exercise 4.9.8 of Kling et al. (2021) Mass Identifications:
# Statistical Methods in Forensic Genetics" Read data and change sex, names to
# agree with previous dataExercise498

KETPex498 = familias2dvir("data-raw/KETPex498.fam", missingPrefix ="MP")

am = KETPex498$am[[1]] |> 
  relabel(new = c("MP1", "MP2", "MP3"),
          old = c("MP3", "MP1", "MP2")) |> 
  swapSex(c("MP2", "MP3")) |> 
  reorderPed(neworder = c(6,5,2,1,3,4))


# Collect and save
KETPex498 = dviData(pm = KETPex498$pm, 
                    am = am, 
                    missing = KETPex498$missing)
usethis::use_data(KETPex498, overwrite = TRUE)

# Check
if(FALSE){
  plotDVI(KETPex498)
  sequentialDVI(KETPex498)$details
  jointDVI(KETPex498)
}


