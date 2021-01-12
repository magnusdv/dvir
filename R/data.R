#' Allele frequencies
#'
#' A dataset containing allele frequencies for 198 markers.
#'
#' @format A list of length 198. Each element is a vector of allele frequencies
#'   for a single marker. The first 27 markers are STRs, the remaining are SNPs.
#'   
"freqsBlind"



#' example1
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

#' example2
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

#' planecrash
#'
#' DVI dataset based on Exercise 3.3 EKM (2015).
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

#' icmp
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
