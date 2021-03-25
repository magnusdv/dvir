#' Data used in the book Kling et al. (2021)
#'
#' Data used in Example 4.8.3 in Kling et al. (2021) 
#' "Mass Identifications: Statistical Methods in Forensic Genetics".
#' There victims are: V1 (female), V2 (male), V3 (female), V4 (male) and V5 (female). 
#' There is one reference family with
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
#'   
#' @examples 
#' \donttest{
#' pm = dataExample483$pm
#' am = dataExample483$am
#' missing = dataExample483$missing
#' 
#' # Plot and find joint solution
#' ncomb(3,2,6,6)
#' 
#' plotPedList(list(pm, am), marker = 1:2, hatched = typedMembers, 
#'             col = list(red = missing))
#' res = jointDVI(pm, am, missing)
#' 
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
#' There are 10 males (V1, V3, ..., V19) and 10 female victims (V2, V4, ..., V20).
#' There are 10 reference families. In each family there is a genotyped grandmother
#' and a missing grandson and a missing granddaughter. The data is simulated according to
#' Vi = Mi, i = 1, ..., 20. 
#' 
#' @format A list of 3 elements:
#'
#'   * `pm`: A list of 20 singletons (victims).
#'
#'   * `am`: A list of 10 pedigrees. 
#'
#'   * `missing`: A vector containing the names of the 20 missing persons.
#'   
#' @examples
#' 
#' \donttest{
#' # Remove comments to run example
#' 
#' pm = sibPairs$pm
#' am = sibPairs$am
#' missing = sibPairs$missing
#' sequentialDVI(pm, am, missing, updateLR = TRUE)
#' # jointDVI(pm, am, missing, threshold = 100)
#' 
#' # Reduce to 15 markers. `sequentialDVI` still gives correct solutions,
#' # but `jointDVI` struggles. Recommend sequential approach or possible to modify the joint?
#' set1 = c("CSF1PO", "D2S1338", "D3S1358", "D5S818", "D7S820", "D8S1179", "D13S317", "D16S539",
#'        "D18S51", "D19S433", "D21S11", "FGA", "TH01", "TPOX", "VWA")
#' 
#' # jointDVI(pm, am, missing, markers = set1, threshold = 10)
#' }
"sibPairs"

