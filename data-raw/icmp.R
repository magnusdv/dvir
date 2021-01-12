load("C:/Users/theg/Dropbox/manuscripts/dvi/dviPaperExamples/icmp.RData")
pm = icmp$pm
am = icmp$am
missing = icmp$MPs
victims = icmp$victims
references = icmp$references
old = labels(am)
new = old
new[c(6, 11, 12, 13, 14, 15, 18:23)] = c("M1","M5", "M2", "M3","M6", "M7",
                                         "M4", "M9", "M8", "M10", "M11", "M12")
am = relabel(am, new, old)
missing = paste0("M", 1:12)
# Check
library(pedtools)
par(mfcol = c(1,1))
plot(am, marker = NULL, col = list(red = missing, blue =references), aff = c(missing, references), 
     deceased = missing)
plotPedList(pm)

# Collect and save
icmp = list(pm = pm, am = am, missing = missing)

usethis::use_data(icmp, overwrite = TRUE)
