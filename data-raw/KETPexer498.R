# Data for Exercise 4.9.8 of "Mass Identifications" (KETP, 2021)

x = familias2dvir("data-raw/KETPex498.fam", missingPrefix = "M")

plotDVI(x)

# Relabel and reorder to agree with book (p.126)
am = x$am[[1]] |>
  relabel(new = c("M1", "M3", "FA1", "FA2"),
          old = c("M3", "M1", "EM1", "added_1")) |> 
  parentsBeforeChildren()

# Collect and save
KETPexer498 = dviData(pm = x$pm, am = am, x$missing)
usethis::use_data(KETPexer498, overwrite = TRUE)

# Check
if(FALSE){
  plotDVI(KETPexer498)
  dviSolve(KETPexer498, detailedOutput = TRUE)
  sequentialDVI(KETPexer498)$details
  dviJoint(KETPexer498)
}


