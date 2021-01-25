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

#' Data. DVI, icmp example
#'
#' DVI dataset based loosely  on 
#' http://www.few.vu.nl/~ksn560/Block-III-PartI-KS-ISFG2017.pdf (page 18).
#' There are 3 female victims, 2 male victims and 6 missing persons of both sexes. 
#' We have  renamed the individuals and simulated data for 13 codis markers.
#'
#' @format A list of 3 elements:
#'
#'   * `pm`: A list of 5 singletons (victims).
#'
#'   * `am`: A pedigree with 12 missing persons
#'
#'   * `missing`: A vector containing the names of the missing persons.
#'   
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
#'   
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
#' pm = sibPairs$pm
#' am = sibPairs$am
#' missing = sibPairs$missing
#' # All correctly identified using all 35 markers:
#' sequentialDVI(pm, am, missing, updateLR = FALSE)
#' 
#' # Six identifications missed if only 15 standard markers are used
#' set1 = c("CSF1PO", "D2S1338", "D3S1358", "D5S818", "D7S820", "D8S1179", "D13S317", "D16S539",
#'        "D18S51", "D19S433", "D21S11", "FGA", "TH01", "TPOX", "VWA")
#' pm15 = selectMarkers(pm, markers = set1)
#' am15 = selectMarkers(am, markers = set1)        
#' sequentialDVI(pm15, am15, missing, updateLR = FALSE)
#' 
#' # All correctly identified using 23 markers:
#' set2 = c("D1S1656", "D2S441", "D10S1248", "D12S391", "D22S1045", "PENTA_D", "PENTA_E", "SE33")
#' pm23 = selectMarkers(pm, markers = c(set1, set2))
#' am23 = selectMarkers(am, markers = c(set1, set2))        
#' sequentialDVI(pm23, am23, missing, updateLR = FALSE)
#' 
#' # I did not wait for below to possibly finish since
#' # 'Assignments to consider in the joint analysis: 2390116
#' res3 = jointDVI(pm, am, missing, verbose = TRUE)
#' 
#' # Works nicely with reduced threshold
#' res3 = jointDVI(pm, am, missing, verbose = TRUE, threshold = 100)
#' 
#' \dontrun{
#' # Error:
#' sequentialDVI(pm, am, missing, updateLR = TRUE, check = FALSE)
#' }
#' 
"sibPairs"

