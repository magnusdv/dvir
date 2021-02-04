# ICMP example: 
# Inspired by page 18 of http://few.vu.nl/~ksn560/Block-III-PartI-KS-ISFG2017.pdf

library(pedtools)

# Reference pedigree
df = read.table(header = TRUE, text = 
"id fid mid sex
  3   0   0   1
  4   0   0   2
 M2   3   4   2
 R2   3   4   1
 R3   3   4   1
 R4   3   4   1
  1   0   0   1
  2   0   0   2
  5   1   2   1
 M1   0   0   2
 R1   1   2   1
 M3   5  M1   1
 M4   5  M1   1
 R5   0   0   2
 M5   5  M1   1
 M6   0   0   2
 M7  M3  M2   2
 R6  M4  R5   1
 M8  M4  R5   1
 M9  M4  R5   2
M10  M5  M6   1
M11  M5  M6   1
M12  M5  M6   2")

refs = c(paste0("R", 1:6))
missing = paste0("M", 1:12)

CODIS = readFreqDatabase("data-raw/codis.txt")

# Simulate PM and AM data
refped = setMarkers(as.ped(df), locusAttributes = CODIS)

am = forrel::profileSim(refped, N = 1, ids = c(refs, missing), seed = 42)[[1]]

pm = list(singleton(id = "V1", sex = 2),
          singleton(id = "V2", sex = 1),
          singleton(id = "V3", sex = 2),
          singleton(id = "V4", sex = 1),
          singleton(id = "V5", sex = 2))
vics = unlist(labels(pm))

# Transfer from true solution
pm = transferMarkers(from = am, to = pm, idsFrom = c("M6", "M10", "M12", "M8", "M1"), idsTo = vics)

# Delete alleles from missing persons in AM
am = setAlleles(am, ids = missing, alleles = 0)

# Collect and save
icmp = list(pm = pm, am = am, missing = missing)

usethis::use_data(icmp, overwrite = TRUE)


#--------------------
# Checks
plot(am, hatched = typedMembers, col = list(red = missing, blue = refs))
dvir::jointDVI(pm, am, missing, verbose = TRUE, numCores = 4)
