#' Plot a DVI problem
#'
#' @param dvi A `dviData` object, typically created with [dviData()].
#' @param pm Either a logical indicating if the PM data should be plotted (as a
#'   set of singletons), or a vector of indices selecting a subset of the PM
#'   samples. Default: TRUE.
#' @param am Either a logical indicating if the AM families data should be
#'   plotted, or a vector of indices selecting a subset of the families.
#'   Default: TRUE.
#' @param hatched A character vector of ID labels, or the name of a function. By
#'   default, typed individuals are hatched.
#' @param col A list of colour vectors (see [pedtools::plot.ped()]). By default,
#'   missing members of `dvi$am` are shown in red.
#' @param frames A logical, by default TRUE.
#' @param titles A character of length 2.
#' @param widths A numeric with relative plot widths.
#' @param nrowPM The number of rows in the array of PM singletons.
#' @param dev.height,dev.width Plot height and widths in inches. These are
#'   optional, and only relevant if `newdev = TRUE`.
#' @param newdev A logical indicating if a new plot window should be opened.
#' @param ... Further parameters to be passed on to [pedtools::plot.ped()],
#'   e.g., `marker`, `cex`, `cex.main`, `symbolsize`.
#'
#' @return NULL
#'
#' @examples
#'
#' plotDVI(example2)
#'
#' # Override default layout of PM singletons
#' plotDVI(example2, nrowPM = 3)
#'
#' # Subset
#' plotDVI(example2, pm = 1:2, am = 1, titles = c("PM (1-2)", "AM (1)"))
#'
#' # AM only
#' plotDVI(example2, pm = FALSE, titles = "AM families")
#'
#' # Further options
#' # plotDVI(example2, new = T, frames = FALSE, marker = 1, cex = 1.2, nrowPM = 1)
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics grconvertX grconvertY layout mtext par rect
#' @export
plotDVI = function(dvi, pm = TRUE, am = TRUE, hatched = typedMembers, col = list(red = dvi$missing), 
                   frames = TRUE, titles = c("PM", "AM"), widths = NULL, nrowPM = NA,
                   dev.height = NULL, dev.width = NULL, newdev = !is.null(c(dev.height, dev.width)), ...) {
  
  # Ensure proper dviData object
  dvi = consolidateDVI(dvi)
  
  PM = if(isTRUE(pm)) dvi$pm else if(!is.logical(pm)) dvi$pm[pm] else NULL
  AM = if(isTRUE(am)) dvi$am else if(!is.logical(am)) dvi$am[am] else NULL
  
  npm = length(PM)
  nam = length(AM)
  if(npm + nam == 0) {
    message("Nothing to plot")
    return(invisible())
  }
    
  if(npm == 0) {
    if(nam == 1)
      plot(AM[[1]], hatched = hatched, col = col, title = titles[length(titles)], ...)
    else
      plotPedList(AM, hatched = hatched, col = col, widths = widths, frames = FALSE, titles = titles[length(titles)],
                  groups = list(1:nam), dev.height = dev.height, dev.width = dev.width, newdev = newdev, ...)
    return(invisible())
  }
  
  if(nam == 0) {
    plotPM(PM, hatched = hatched, col = col, title = titles[1], nrow = nrowPM, ...)
    return(invisible())
  }
  
  # Calculate widths
  alignlist = lapply(AM, .pedAlignment)
  maxGen = max(vapply(alignlist, function(a) a$maxlev, 1))
  pmDims = pmArrayDim(length(PM), nrow = nrowPM)
  nrpm = pmDims[1]
  ncpm = pmDims[2]
  
  if(is.null(widths)) {
    xranges = vapply(alignlist, function(a) max(1, diff(a$xrange)), 1)
    widths = sqrt(xranges - 1) + 1
    pmwid = max(sqrt(ncpm), sum(widths)/3)
    widths = c(pmwid, widths)
  }
  
  # Layout of plot regions
  if (newdev) {
    dev.height = dev.height %||% {max(3, 1.2 * maxGen) + 0.3}
    dev.width = dev.width %||% (sum(widths + 1))
    dev.new(height = dev.height, width = dev.width, noRStudioGD = TRUE)
  }
  
  # Set graphical parameters (include mfrow to ensure layout is reverted on exit)
  opar = par(oma = c(0, 0, 3, 0), xpd = NA, mfrow = c(1,1), mar = c(0,0,0,0))
  on.exit(par(opar))
  
  layout(rbind(1:(1 + length(AM))), widths = widths)
  
  plotPM(PM, nrow = nrpm, hatched = hatched, col = col, margins = 2, ...)
  for(a in AM)
    plot(a, hatched = hatched, col = col, margins = 2, ...)
  
  # X coordinate of each plot region (converted to value in [0,1]).
  ratios = c(0, widths[1], sum(widths))/sum(widths)
  
  # Draw frames
  if(frames) {
    
    margIn = grconvertY(1, from = "lines", to = "inches")
    margx = grconvertX(margIn, from = "inches", to = "ndc")
    margy = grconvertY(margIn, from = "inches", to = "ndc")
    
    rect(xleft   = grconvertX(ratios[1:2] + margx, from = "ndc"),
         ybottom = grconvertY(1 - margy, from = "ndc"),
         xright  = grconvertX(ratios[2:3] - margx, from = "ndc"),
         ytop    = grconvertY(margy, from = "ndc"), 
         xpd = NA)
  }
  
  # Add titles
  midpoints = ratios[1:2] + diff(ratios)/2
  args = list(...)
  cex.tit = cex = args[["cex.main"]] %||% args[["cex"]]
  mtext(titles, outer = TRUE, at = midpoints, cex = cex.tit, font = 2)
}


pmArrayDim = function(N, nrow = NA) {
  if(is.na(nrow))
    nr = if(N <= 3) N else ceiling(sqrt(N))
  else 
    nr = min(N, nrow)
  c(nr, ceiling(N/nr))
}


plotPM = function(pm, nrow = NA, ...) {

  # alignment ---------------------------------------------------------------
  
  N = length(pm)
  dims = pmArrayDim(N, nrow)
  nr = dims[1]
  nc = dims[2]
  
  nid = matrix(c(1:N, rep(0, nr*nc - N)), nrow = nr, ncol = nc, byrow = TRUE)
  pos = col(nid) - 1
  
  n = rowSums(nid > 0)
  fam = spouse = matrix(0, nrow = nr, ncol = nc)
  
  plist = list(n = n, nid = nid, pos = pos, fam = fam, spouse = spouse)
  align = pedtools::.pedAlignment(pm, plist = plist)
  
  # annotation --------------------------------------------------------------
  
  annotlist = lapply(pm, function(s) .pedAnnotation(s, ...))
  annot1 = annotlist[[1]]
  
  mrg = function(key, def)
    unlist(lapply(annotlist, function(a) a[[key]] %||% def))
  
  annot = list(
    title = annot1$title,
    textUnder = mrg("textUnder", ""),
    textAbove = mrg("textAbove", ""),
    textInside = mrg("textInside", ""),
    colvec = mrg("colvec", 1),
    aff01 = mrg("aff01", 0),
    density = annot1$density,
    angle = annot1$angle,
    carrierTF = mrg("carrierTF", FALSE),
           deceasedTF = mrg("deceasedTF", FALSE))
  
  p = drawPed(align, annot, ...)
  invisible(p)
}



#' Plot DVI solution
#'
#' @param dvi A dviData object
#' @param assignment A named character of the format `c(victim = missing, ...)`,
#'   or a data frame produced by [jointDVI()].
#' @param k An integer; the row number when `assignment` is a data frame.
#' @param ... Parameters passed on to [plotDVI()].
#'
#' @return NULL.
#'
#' @examples
#'
#' res = jointDVI(example2, verbose = FALSE)
#'
#' plotSolution(example2, res)
#' plotSolution(example2, res, marker = 1)
#'
#' # Non-optimal solutions
#' plotSolution(example2, res, k = 2, pm = FALSE)
#' plotSolution(example2, res, k = 2, cex = 1.3)
#'
#' @export
plotSolution = function(dvi, assignment, k = 1, ...) {
  
  # Ensure proper dviData object
  dvi = consolidateDVI(dvi)
  
  a = assignment
  
  if(is.data.frame(a))
    a = unlist(a[k, !names(a) %in% c("loglik", "LR", "posterior"), drop = FALSE])
  
  mtch = a[a != "*"]
  newlabs = paste(names(mtch), mtch, sep = "=")
  refs = typedMembers(dvi$am)
  
  if(length(mtch)) {
    am = relabel(dvi$am, old = mtch, new = newlabs)
    am = transferMarkers(dvi$pm[names(mtch)], am, idsFrom = names(mtch), idsTo = newlabs, erase = FALSE)
    dvi$am = am
  }
  
  dvi$pm = dvi$pm[a == "*"]
  
  stillmissing = setdiff(dvi$missing, mtch)
  
  plotDVI(dvi, fill = list("green" = newlabs), hatched = refs, col = list("red" = stillmissing), 
          lwd = list("2" = c(stillmissing)), ...)
}