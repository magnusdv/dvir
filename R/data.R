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
#'   * `MPs`: A vector containing the names of the missing persons.
#'   
"example1"


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
#'   * `MPs`: A vector containing the names of the missing persons.
#'   
"planecrash"
