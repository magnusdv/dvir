#' Plot a DVI problem
#'
#' @param dvi A `dviData` object, typically created with [dviData()].
#' @param hatched A character vector of ID labels, or the name of a function. By
#'   default, typed individuals are hatched.
#' @param col A list of colour vectors (see [pedtools::plot.ped()]). By default,
#'   missing members of `dvi$am` are shown in red.
#' @param frames A logical, by default TRUE.
#' @param titles A character of length 2.
#' @param widths A numeric with relative plot widths.
#' @param dev.height,dev.width Plot height and widths in inches. These are
#'   optional, and only relevant if `newdev = TRUE`.
#' @param newdev A logical indicating if a new plot window should be opened.
#' @param ... Further parameters to be passed on to [pedtools::plot.ped()], e.g.,
#'   `marker`, `cex`, `cex.main`, `symbolsize`.
#'
#' @return NULL
#'
#' @examples
#'
#' plotDVI(example1, marker = 1)
#'
#' # plotDVI(icmp)
#'
#' # plotDVI(example2, new = T, frames = FALSE, marker = 1, cex = 1.2)
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics grconvertX grconvertY layout mtext par rect
#' @export
plotDVI = function(dvi, hatched = typedMembers, col = list(red = dvi$missing), 
                   frames = TRUE, titles = c("PM", "AM"), widths = NULL,
                   dev.height = NULL, dev.width = NULL,
                   newdev = !is.null(c(dev.height, dev.width)), ...) {
  
  dvi = consolidate(dvi)
  
  pm = dvi$pm
  am = dvi$am
  
  # Calculate widths
  alignlist = lapply(am, .pedAlignment)
  maxGen = max(vapply(alignlist, function(a) a$maxlev, 1))
  nrpm = maxGen
  ncpm = ceiling(length(pm)/nrpm)
  
  if(is.null(widths)) {
    xranges = vapply(alignlist, function(a) max(1, diff(a$xrange)), 1)
    widths = sqrt(xranges - 1) + 1
    pmwid = max(sqrt(ncpm), sum(widths)/3)
    widths = c(pmwid, widths)
  }
  
  # Layout of plot regions
  if (newdev) {
    dev.height = dev.height %||% {max(3, 1 * maxGen) + 0.3}
    dev.width = dev.width %||% (sum(widths + 1))
    dev.new(height = dev.height, width = dev.width, noRStudioGD = TRUE)
  }
  
  # Set graphical parameters (include mfrow to ensure layout is reverted on exit)
  opar = par(oma = c(0, 0, 3, 0), xpd = NA, mfrow = c(1,1), mar = c(0,0,0,0))
  on.exit(par(opar))
  
  layout(rbind(1:(1 + length(am))), widths = widths)
  
  plotPM(pm, nrow = nrpm, hatched = hatched, col = col, margins = 2, ...)
  for(a in am)
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


plotPM = function(pm, nrow = NA, ...) {

  # alignment ---------------------------------------------------------------
  
  N = length(pm)
  nr = if(is.na(nrow)) floor(sqrt(N)) else min(N, nrow)
  nc = ceiling(N/nr)
  
  nid = matrix(c(1:N, rep(0, nr*nc - N)), nrow = nr, ncol = nc)
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
  
  annot = list(textUnder = mrg("textUnder", ""),
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
