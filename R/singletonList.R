#' Creates a list of singletons
#' 
#' Wrapper for `singleton`` function of pedtools. Creates a list of singletons 
#' with default individual and family ids and alleles sampled from a SNP marker with equally frequent alleles
#' 1 and 2.
#' 
#' @param n integer Noof singletons
#' @param ids Character vector. Names of victims
#' @param famid Character vector. Names of families
#' @param als Alleles
#' @param freq Allele frequencies
#' @param sex Integer vector
#' @param amat Matrix A nx2 matrix with alleles sampled from and 2
#' @param seed Integer
#' @return  A list of singletons
#' @examples 
#' singletonList()
#' @export
singletonList = function(n = 2, ids = paste("V", 1:n, sep = ""), 
                         famid = paste("F", 1:n, sep = ""), sex = rep(1,n),
                         als = 1:2, freq = c(0.5, 0.5),
                         amat = matrix(sample(als, 2*n, rep = TRUE), ncol = 2),
                         seed = NULL){
  set.seed(seed)
  singles = list()
  for (i in 1:n){
    singles[[i]] = singleton(ids[i], sex = sex[i], famid = famid[i])
    m = marker(singles[[i]], alleles = als, afreq = freq)
    singles[[i]] = addMarkers(singles[[i]], m)
  }
  singles = setAlleles(singles, markers = 1, alleles = amat)
  singles
}
