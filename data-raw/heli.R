## code to prepare `heli` dataset goes here

heli = readRDS("data-raw/heli.rds")
heli
plotDVI(heli)

usethis::use_data(heli, overwrite = TRUE)
