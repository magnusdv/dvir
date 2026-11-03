#' LR power for DVI reference families
#'
#' Simulates the pairwise LR for a selected missing person versus an unrelated candidate.
#' All other missing persons remain untyped in the AM pedigree. The two LR hypotheses
#' always retain the relationships among AM references. PM data, if present in `dvi`, are
#' ignored.
#'
#' @param dvi A `dviData` object. PM profiles are not required.
#' @param fams AM families to include, by name. Default: all.
#' @param missing Missing persons to include. Default: all in the selected families.
#' @param conditional If TRUE (default), keep reference genotypes fixed. If FALSE,
#'   simulate them with their original missing-data pattern.
#' @param simulate Generating hypothesis: "both" (default), "H1" (true missing person), or
#'   "H2" (unrelated person).
#' @param markers Markers to use. Default: all attached AM markers.
#' @param nsim Number of simulations per hypothesis and missing person.
#' @param threshold Positive LR threshold for the summary statistics.
#' @param seed Optional RNG seed; applied once for the entire calculation.
#' @param verbose A logical.
#'
#' @return A list with `results` (one named H1/H2 `LRpowerResult` list per family–MP
#'   pairing), `summary` (one row per pairing), `labels` and `params`. Individual results
#'   may be plotted with `forrel::LRpowerPlot(data = ...)`.
#'
#' @examples
#' # Simulations use AM pedigrees only; no PM profiles are needed.
#' d = dviData(pm = NULL, am = example2$am, missing = example2$missing)
#' z = dviPower(d, fams = "F2", nsim = 10, seed = 123, verbose = FALSE)
#' forrel::LRpowerPlot(data = z$results[[1]], labels = c("True MP", "Unrelated"))
#'
#' @export
#' @importFrom pedtools relabel singleton getSex removeGenotypes
#' @importFrom forrel LRpower

dviPower = function(dvi, fams = NULL, missing = NULL, conditional = TRUE,
                     simulate = c("both", "H1", "H2"), markers = NULL,
                     nsim = 500, threshold = 1e4, seed = NULL, verbose = TRUE) {

  if(!inherits(dvi, "dviData"))
    stop2("`dvi` must be a dviData object")
  dvi = consolidateDVI(dvi)
  simulate = match.arg(simulate)
  hyp = if(simulate == "both") c("H1", "H2") else simulate
  fams = fams %||% names(dvi$am)
  bad = setdiff(fams, names(dvi$am))
  if(length(bad)) stop2("Unknown AM family: ", bad)
  missing = unique(as.character(missing %||% dvi$missing))
  bad = setdiff(missing, dvi$missing)
  if(length(bad)) stop2("Unknown missing person: ", bad)

  # One calculation per selected MP, including multiple MPs in complex families.
  family = getFamily(dvi, ids = missing)
  keep = family %in% fams
  pairs = data.frame(Family = unname(family[keep]), Missing = missing[keep])
  if(!nrow(pairs))
    stop2("No missing persons in the selected families")
  if(!is.null(threshold) && (length(threshold) != 1L || !is.finite(threshold) || threshold <= 0))
    stop2("`threshold` must be a positive number or NULL")

  results = lapply(seq_len(nrow(pairs)), function(i) {
    fam = pairs$Family[i]
    mp = pairs$Missing[i]
    ped = dvi$am[[fam]]
    if(!nMarkers(ped))
      stop2("No markers in family ", fam)
    ped = removeGenotypes(ped, ids = intersect(dvi$missing, ped$ID))
    poi = make.unique(c(ped$ID, "_POI_"))[length(ped$ID) + 1L]

    # H1: candidate is MP. H2: candidate unrelated, references unchanged.
    H1 = relabel(ped, old = mp, new = poi)
    H2 = list(ped, singleton(poi, sex = getSex(ped, mp)))
    z = lapply(seq_along(hyp), function(j) {
      h = hyp[j]
      if(verbose) message(sprintf("%s / %s: %s", fam, mp, h))
      forrel::LRpower(H1, H2, truePed = if(h == "H1") H1 else H2,
                      ids = poi, markers = markers, source = "numerator",
                      nsim = nsim, threshold = threshold, conditional = conditional,
                      seed = if(i == 1L && j == 1L) seed else NULL,
                      verbose = verbose)
    })
    names(z) = hyp
    z
  })
  names(results) = paste(pairs$Family, pairs$Missing, sep = ":")

  # Use log10 LRs for summaries; no resimulation is needed to change the threshold.
  probability = function(z, h) {
    if(is.null(z[[h]]) || is.null(threshold)) return(NA_real_)
    mean(z[[h]]$log10LRperSim >= log10(threshold))
  }
  med = function(z, h) {
    if(is.null(z[[h]])) return(NA_real_)
    stats::median(z[[h]]$log10LRperSim)
  }
  pairs$P_H1 = vapply(results, probability, numeric(1), h = "H1")
  pairs$P_H2 = vapply(results, probability, numeric(1), h = "H2")
  pairs$MedianLog10LR_H1 = vapply(results, med, numeric(1), h = "H1")
  pairs$MedianLog10LR_H2 = vapply(results, med, numeric(1), h = "H2")

  list(results = results, summary = pairs,
       labels = c("True missing person", "Unrelated person"),
       params = list(conditional = conditional, simulate = simulate, markers = markers,
                     nsim = nsim, threshold = threshold, seed = seed))
}
