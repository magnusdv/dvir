#' DVI dataset: Generational trio
#'
#' A proof-of-concept dataset involving three missing members (child, father,
#' grandfather) of a single family. With the given data, stepwise victim
#' identification fails to find the correct solution, while joint identification
#' succeeds.
#'
#' @format A `dviData` object with the following content:
#'
#'   * `pm`: A list of 3 singletons (victims).
#'
#'   * `am`: A pedigree with three missing persons and one typed reference
#'   individual.
#'
#'   * `missing`: A vector containing the names of the missing persons.
#'
#' @examples    
#' example1
#' 
#' plotDVI(example1, marker = 1)
#' 
#' jointDVI(example1)  
#' 
"example1"

#' DVI dataset: Two reference families
#'
#' A small DVI example with three victims, and three missing persons from two reference
#' families
#'
#' @format A `dviData` object with the following content:
#'
#'   * `pm`: A list of 3 singletons (victims).
#'
#'   * `am`: A list of 2 pedigrees with three missing persons and one typed reference
#'   individual.
#'
#'   * `missing`: A vector containing the names of the missing persons.
#'   
#' @examples    
#' example2
#' 
#' plotDVI(example2, marker = 1, nrowPM = 3)
#' 
#' jointDVI(example2)
#' 
"example2"

#' DVI dataset: Simulated plane crash
#'
#' A simulated dataset based on Exercise 3.3 in Egeland et al. "Relationship
#' Inference with Familias and R" (2015).
#'
#' The 15 markers are `CSF1PO`, `D13S317`, `D16S539`, `D18S51`, `D21S11`,
#' `D3S1358`, `D5S818`, `D7S820`, `D8S1179`, `FGA`, `PENTA_D`, `PENTA_E`,
#' `TH01`,  `TPOX`,   and  `VWA`.
#'
#' Source code for the simulation, and a file containing the allele frequencies,
#' can be found in the `data-raw` folder of the GitHub repository:
#' https://github.com/magnusdv/dvir.
#'
#' @format A `dviData` object with the following content:
#'
#'   * `pm`: A list of 8 female singletons (victims).
#'
#'   * `am`: A list of 5 pedigrees, each with one missing member and one
#'   genotyped member.
#'
#'   * `missing`: A vector containing the names of the missing persons.
#'
#' @examples
#' planecrash
#' 
#' # plotDVI(planecrash)
#' 
#' # Markers and allele frequencies
#' db = pedtools::getFreqDatabase(planecrash$pm)
#' db 
#' 
"planecrash"


#' DVI dataset: A large reference pedigree
#'
#' DVI dataset based loosely on the ICMP workshop material
#' https://www.few.vu.nl/~ksn560/Block-III-PartI-KS-ISFG2017.pdf (page 18).
#' There are 3 female victims, 2 male victims and 6 missing persons of both
#' sexes. We have renamed the individuals and simulated data for 13 CODIS
#' markers (see Details).
#'
#' The 13 markers are, in order: `CSF1PO`, `D3S1358`, `D5S818`,`D7S820`,
#' `D8S1179`, `D13S317`, `D16S539`, `D18S51`, `D21S11`, `FGA`, `TH01`, `TPOX`,
#' and `vWA`.
#'
#' Source code for the simulation, and a file containing the allele frequencies,
#' can be found in the `data-raw` folder of the GitHub repository:
#' https://github.com/magnusdv/dvir.
#'
#' @format A `dviData` object with the following content:
#'
#'   * `pm`: A list of 5 singletons (victims).
#'
#'   * `am`: A reference pedigree with 6 genotyped members and 12 missing
#'   persons.
#'
#'   * `missing`: A vector containing the names of the missing persons.
#'
#'
#' @examples
#' icmp
#' 
#' # plotDVI(icmp)
#' 
#' # Markers and allele frequencies
#' db = pedtools::getFreqDatabase(icmp$pm)
#' db
#' 
"icmp"


#' DVI dataset: Family grave
#'
#' Family grave data in Kling et al. (2021) "Mass Identifications: Statistical
#' Methods in Forensic Genetics". There are 5 female victims and 3 male victims.
#' There is one reference family with 5 missing females and 3 missing males.
#' There are 23 markers, no mutation model.
#'
#'
#' @format A `dviData` object with the following content:
#'
#'   * `pm`: A list of 8 singletons (victims).
#'
#'   * `am`: A pedigree with 8 missing persons.
#'
#'   * `missing`: A vector containing the names of the missing persons.
#'
#' @examples
#' grave
#' 
#' # plotDVI(grave, marker = 1)
#'
#' # jointDVI(grave)
#' 
"grave"


#' Data used in the book Kling et al. (2021)
#'
#' Data used in last example of Chapter 4 in Kling et al. (2021) "Mass
#' Identifications: Statistical Methods in Forensic Genetics". There are 2
#' female victims, 2 male victims. There are four reference families with 2
#' missing females and 2 missing males. There are 21 markers. An `equal mutation
#' mode with rate 0.005 is specified.
#'
#'
#' @format A `dviData` object with the following content:
#'
#'   * `pm`: A list of 4 singletons (victims).
#'
#'   * `am`: A list of 3 pedigrees.
#'
#'   * `missing`: A vector containing the names of the missing persons.
#'
#' @examples
#' KETPch4
#' 
#' plotDVI(KETPch4, nrowPM = 4)
#' 
#' 
"KETPch4"


#' Data used in the book Kling et al. (2021)
#'
#' Data used in Exercise 4.9.7 in Kling et al. (2021) "Mass Identifications:
#' Statistical Methods in Forensic Genetics". There are 3 female victims and 3
#' reference families with 3 missing females. There are 23 markers, equal
#' mutation model, rate 0.001.
#'
#'
#' @format A `dviData` object with the following content:
#'
#'   * `pm`: A list of 3 singletons (victims).
#'
#'   * `am`: A list of 3 pedigrees.
#'
#'   * `missing`: A vector containing the names of the missing persons.
#'
#' @examples
#' plotDVI(KETPex497, nrowPM = 3)
#'
#'
"KETPex497"


#' Data used in the book Kling et al. (2021)
#'
#' Data used in Exercise 4.9.8 in Kling et al. (2021) "Mass Identifications:
#' Statistical Methods in Forensic Genetics". There are 2 female victims and 1
#' male. There is one reference family with 2 missing females and one missing
#' male. There are 16 markers, equal mutation model, rate 0.001.
#'
#'
#' @format A `dviData` object with the following content:
#'
#'   * `pm`: A list of 3 singletons (victims).
#'
#'   * `am`: A list of 1 pedigree.
#'
#'   * `missing`: A vector containing the names of the missing persons.
#'   
#' @examples 
#' 
#' plotDVI(KETPex498, nrowPM = 3)
#' 
"KETPex498"


#' Data used in the book Kling et al. (2021)
#'
#' Data used in Example 4.8.1 in Kling et al. (2021) "Mass Identifications:
#' Statistical Methods in Forensic Genetics". There  victims are V1 and V2, both
#' females. There is one reference family with 2 missing persons, both females.
#' There are 21 markers, no mutation model.
#'
#' @format A `dviData` object with the following content:
#'
#'   * `pm`: A list of 2 singletons (victims).
#'
#'   * `am`: A list of 1 pedigree.
#'
#'   * `missing`: A vector containing the names of the missing persons.
#'
#' @examples
#' plotDVI(KETPex481, marker = 1)
#' 
"KETPex481"


#' Dataset: Exclusion example
#'
#' This data is based on a real case, but pedigrees have been changed and
#' marker data simulated to preserve anonymity.
#' 
#' @format A `dviData` object with the following content:
#'
#'   * `pm`: A list of 16 singletons (male victims).
#'
#'   * `am`: A list of 15 pedigrees, each with one missing person 
#'
#'   * `missing`: A vector containing the names of the 15 missing persons.
#'   
#' @examples
#' 
#' exclusionExample
#' 
"exclusionExample"
