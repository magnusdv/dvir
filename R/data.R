#' Allele frequencies
#'
#' A dataset containing allele frequencies for 198 markers.
#'
#' @format A list of length 198. Each element is a vector of allele frequencies
#'   for a single marker. The first 27 markers are STRs, the remaining are SNPs.
#'   
"freqsBlind"



#' Data. Stepwise identification may fail
#'
#' A proof-of-concept dataset showing that showing that stepwise victim
#' identification may fail
#'
#' @format A list of 3 elements:
#'
#'   * `pm`: A list of 3 singletons (victims).
#'
#'   * `am`: A pedigree with three missing persons and one typed reference
#'   individual.
#'
#'   * `missing`: A vector containing the names of the missing persons.
#'   
"example1"

#' Data. DVI with two reference families
#'
#' A small DVI example with three victims, and three missing persons from two reference
#' families
#'
#' @format A list of 3 elements:
#'
#'   * `pm`: A list of 3 singletons (victims).
#'
#'   * `am`: A list of 2 pedigrees with three missing persons and one typed reference
#'   individual.
#'
#'   * `missing`: A vector containing the names of the missing persons.
#'   
#' @examples    
#' pm = example2$pm
#' am = example2$am
#' missing = example2$missing
#' jointDVI(pm, am, missing,)
"example2"

#' Data. planecrash
#'
#' DVI dataset based on simulated data in Exercise 3.3 Egeland et al. (2015).
#'
#' @format A list of 3 elements:
#'
#'   * `pm`: A list of 8 female singletons (victims).
#'
#'   * `am`: A list of 5 pedigrees, each with one missing member and one
#'   genotyped member.
#'
#'   * `missing`: A vector containing the names of the missing persons.
#'   
"planecrash"

#' DVI example with large reference pedigree
#'
#' DVI dataset based loosely on the ICMP workshop material
#' http://www.few.vu.nl/~ksn560/Block-III-PartI-KS-ISFG2017.pdf (page 18).
#' There are 3 female victims, 2 male victims and 6 missing persons of both sexes. 
#' We have renamed the individuals and simulated data for 13 codis markers.
#'
#' @format A list of 3 elements:
#'
#'   * `pm`: A list of 5 singletons (victims).
#'
#'   * `am`: A reference pedigree with 6 genotyped members and 12 missing persons
#'
#'   * `missing`: A vector containing the names of the missing persons.
#'   
#' @examples 
#' library(pedtools)
#' 
#' plot(icmp$am, hatched = typedMembers, col = list(red = icmp$missing))
#' 
#' \donttest{
#' jointDVI(icmp$pm, icmp$am, icmp$missing)
#' }   
"icmp"

#' Data. DVI, large pedigree
#'
#' Family grave data in Kling et al. (2021) 
#' "Mass Identifications: Statistical Methods in Forensic Genetics"
#' There are 5 female victims, 3 male victims. There is one
#' reference family with 5 missing females and 3 missing males. 
#' There are 23 markers, no mutation model.
#' 
#'
#' @format A list of 3 elements:
#'
#'   * `pm`: A list of 8 singletons (victims).
#'
#'   * `am`: A pedigree with 8 missing persons
#'
#'   * `missing`: A vector containing the names of the missing persons.
#'   
#' @examples   
#' pm = grave$pm # The list of missing persons
#' am = grave$am # The reference family pedigree
#' missing = grave$missing # The names of the missing persons
#' library(forrel)
#' plot(am, marker = 1)
#' plotPedList(pm, frames = FALSE, marker = 1)
"grave"

#' Data. Used in the book Kling et al. (2021)
#'
#' Data used in last example of Chapter 4  in Kling et al. (2021) 
#' "Mass Identifications: Statistical Methods in Forensic Genetics"
#' There are 2 female victims, 2 male victims. There are four
#' reference families with 2 missing females and 2 missing males. 
#' There are 21 markers, equal mutation model, rate 0.005.
#' 
#'
#' @format A list of 3 elements:
#'
#'   * `pm`: A list of 4 singletons (victims).
#'
#'   * `am`: A list of 3 pedigrees 
#'
#'   * `missing`: A vector containing the names of the missing persons.
#'  @examples 
#'  \donttest{
#'   pm = dataCh4$pm
#'   am = dataCh4$am
#'   missing = dataCh4$missing
#'   res = jointDVI(pm, am, missing, disableMutations = FALSE)
#'   head(res[c(1, 2, 30, 49),])
#'   }
"dataCh4"

#' Data. Exercise, Kling et al. (2021)
#'
#' Data used in Exercise 4.9.7 in Kling et al. (2021) 
#' "Mass Identifications: Statistical Methods in Forensic Genetics"
#' There are 3 female victims and 3 
#' reference families with 3 missing females. 
#' There are 23 markers, equal mutation model, rate 0.001.
#' 
#'
#' @format A list of 3 elements:
#'
#'   * `pm`: A list of 3 singletons (victims).
#'
#'   * `am`: A list of 3 pedigrees 
#'
#'   * `missing`: A vector containing the names of the missing persons.
#'   
"dataExercise497"

#' Data. Exercise, Kling et al. (2021)
#'
#' Data used in Exercise 4.9.8 in Kling et al. (2021) 
#' "Mass Identifications: Statistical Methods in Forensic Genetics"
#' There are 2 female victims and 1 male, 
#' There is one reference family with 2 missing females and one missing male. 
#' There are 16 markers, equal mutation model, rate 0.001.
#' 
#'
#' @format A list of 3 elements:
#'
#'   * `pm`: A list of 3 singletons (victims).
#'
#'   * `am`: A list of 1 pedigree 
#'
#'   * `missing`: A vector containing the names of the missing persons.
#'   
"dataExercise498"

#' Data. Exercise, Kling et al. (2021)
#'
#' Data used in Example 4.8.1 in Kling et al. (2021) 
#' "Mass Identifications: Statistical Methods in Forensic Genetics"
#' There are victims V1 (female)  V2 (female). 
#' There is one reference family with missing persons
#' 2 missing persons, both females..
#' There are 21 markers,no mutation model.
#' 
#'
#' @format A list of 3 elements:
#'
#'   * `pm`: A list of 2 singletons (victims).
#'
#'   * `am`: A list of 1 pedigree 
#'
#'   * `missing`: A vector containing the names of the missing persons.
#' @examples 
#' \donttest{
#' pm = dataExample481$pm
#' am = dataExample481$am
#' missing = dataExample481$missing
#' ncomb(2, 2, 0, 0)
#' windows()
#' plotPedList(list(pm, am), marker = 1:2, hatched = typedMembers, col = list(red = missing))
#' jointDVI(pm, am, missing, verbose = FALSE)
#' }
"dataExample481"

#' Data. Exercise, Kling et al. (2021)
#'
#' Data used in Example 4.8.3 in Kling et al. (2021) 
#' "Mass Identifications: Statistical Methods in Forensic Genetics"
#' There are victims V1 (female), V2 (male), V3 (female), V4 (male) and V5 (female). 
#' There is one reference family with missing persons
#' 12 missing persons, 6 females and 6 males.
#' There are 13 codis markers, equal mutation model, rate 0.001.
#' 
#'
#' @format A list of 3 elements:
#'
#'   * `pm`: A list of 5 singletons (victims).
#'
#'   * `am`: A list of 1 pedigree 
#'
#'   * `missing`: A vector containing the names of the missing persons.
#' @examples 
#' \donttest{
#' pm = dataExample483$pm
#' am = dataExample483$am
#' missing = dataExample483$missing
#' ncomb(3,2,6,6)
#' # Make large window prior to plotting, e.g. by windows()
#' plotPedList(list(pm, am), marker = 1:2, hatched = typedMembers, 
#'             col = list(red = missing))
#' res = jointDVI(pm, am, missing)
#' # Transfer victims and check solution
#' amId = transferMarkers(pm, am, idsFrom = paste0("V",1:5), 
#'                        idsTo = missing[1:5], erase = FALSE)
#' plot(amId, marker = 1:2, hatched = typedMembers, 
#'             col = list(red = missing))
#' }
"dataExample483"

#' Data. Simulated sib pairs
#'
#' The purpose of this data is to challenge brute force methods. We use the
#' the database NorwegianFrequencies.
#' There are 10 male (V1, V3, ..., V19) and 10 female victims (V2, V4, ..., V20).
#' There are 10 reference families. In each family there is a genotyped grandmother
#' and a missing grandson and a missing granddaughter. The data is simulated according to
#' Vi = Mi, i = 1, ..., 20. 
#' 
#' @format A list of 3 elements:
#'
#'   * `pm`: A list of 20 singletons (victims).
#'
#'   * `am`: A list of 10 pedigrees 
#'
#'   * `missing`: A vector containing the names of the 20 missing persons.
#'   
#' @examples 
#' \donttest{
#' pm = sibPairs$pm
#' am = sibPairs$am
#' missing = sibPairs$missing
#' sequentialDVI(pm, am, missing, updateLR = TRUE)
#' jointDVI(pm, am, missing, threshold = 100)
#' 
#' # Reduce to 15 markers. `sequentialDVI` still gives correct solutions,
#' # but `jointDVI` struggles. Recommend sequential approach or possible to modify the joint?
#' set1 = c("CSF1PO", "D2S1338", "D3S1358", "D5S818", "D7S820", "D8S1179", "D13S317", "D16S539",
#'        "D18S51", "D19S433", "D21S11", "FGA", "TH01", "TPOX", "VWA")
#' 
#' jointDVI(pm, am, missing, markers = set1, threshold = 10)
#' }
"sibPairs"

