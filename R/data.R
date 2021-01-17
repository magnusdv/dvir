#' Allele frequencies
#'
#' A dataset containing allele frequencies for 198 markers.
#'
#' @format A list of length 198. Each element is a vector of allele frequencies
#'   for a single marker. The first 27 markers are STRs, the remaining are SNPs.
#'   
"freqsBlind"



#' Data showing that stepwise identification may fail
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

#' Tutorial DVI data with two reference families
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

#' planecrash data
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

#' icmp data
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

#' grave data
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

#' dataCh4: data used in the book Kling et al. (2021)
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
