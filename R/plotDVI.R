#' Plot a DVI problem
#'
#' @param dvi A `dviData` object, typically created with [dviData()].
#' @param pm Either a logical indicating if the PM data should be plotted (as a
#'   set of singletons), or a vector of indices selecting a subset of the PM
#'   samples. Default: TRUE.
#' @param am Either a logical indicating if the AM families data should be
#'   plotted, or a vector of indices selecting a subset of the families.
#'   Default: TRUE.
#' @param style An integer (currently 1 or 2) indicating the style of the plot.
#' @param hatched A character vector of ID labels, or the name of a function. By
#'   default, typed individuals are hatched.
#' @param fill,cex,col,lwd,carrier Arguments passed on to [pedtools::plot.ped()].
#' @param frames A logical, by default TRUE.
#' @param titles A character of length 2.
#' @param famnames A logical. If NA (default) family names are included if there
#'   are multiple families.
#' @param widths,heights Numeric with relative columns widths / row heights, to
#'   be passed on to `layout()`.
#' @param nrowPM The number of rows in the array of PM singletons.
#' @param nrowAM The number of rows in the array of AM families.
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
#' plotDVI(example2, nrowPM = 1)
#'
#' # Subset
#' plotDVI(example2, pm = 1:2, am = 1, titles = c("PM (1-2)", "AM (1)"))
#'
#' # AM only
#' plotDVI(example2, pm = FALSE, titles = "AM families")
#'
#' # Further plot options
#' plotDVI(example2, frames = FALSE, marker = 1, cex = 1.2, nrowPM = 3, nrowAM = 2,
#'   textAnnot = list(inside = list(c(M1 = "?", M2 = "?", M3 = "?"), cex = 1.5)))
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics grconvertX grconvertY layout mtext par rect
#' @export
plotDVI = function(dvi, pm = TRUE, am = TRUE, style = 1, famnames = NA,  
                   hatched = typedMembers, frames = TRUE, titles = c("PM", "AM"),
                   cex = NA, col = 1, lwd = 1, fill = NA, carrier = NULL, 
                   widths = NULL, heights = NULL, nrowPM = NA, nrowAM = NA, 
                   dev.height = NULL, dev.width = NULL, 
                   newdev = !is.null(c(dev.height, dev.width)), ...) {
  
  # Ensure proper dviData object
  dvi = consolidateDVI(dvi)
  
  # Extra args
  args = list(...)
  
  PM = if(isTRUE(pm)) dvi$pm else if(!is.logical(pm)) dvi$pm[pm] else NULL
  AM = if(isTRUE(am)) dvi$am else if(!is.logical(am)) dvi$am[am] else NULL
  
  npm = length(PM)
  nam = length(AM)
  if(npm + nam == 0) {
    message("Nothing to plot")
    return(invisible())
  }
    
  if(nam == 0) {
    if(newdev)
      stop2("Cannot plot only PM to new device")
    plotPM(PM, hatched = hatched, title = titles[1], nrow = nrowPM, ...)
    return(invisible())
  }
  
  # Ad hoc cex increase for > 2 panels
  if(is.na(cex)) {
    npanels = nam + (npm > 0)
    cex = if(npanels %in% 3:6) 1.3 else if(npanels %in% 7:9) 1.2 else 1
  }
  
  if(is.na(famnames))
    famnames = nam > 1
  
  # AM layout
  amDim = amArrayDim(nam, nrowAM)
  pmDim = pmArrayDim(npm, nrowPM)
  layoutMat = matrix(1:prod(amDim), nrow = amDim[1], ncol = amDim[2], byrow = TRUE)
  
  if(npm > 0)
    layoutMat = cbind(1, layoutMat + 1)
  
  # Relative widths/heights & plot dimensions
  plotdims = findPlotDims(AM, amDim, pmDim)
  widths = widths %||% plotdims$widths
  heights = heights %||% plotdims$heights |> (\(x) {x[1] = x[1] + famnames * 0.2; x})()

  # Open new device?
  if (newdev) 
    dev.new(height = dev.height %||% plotdims$H, 
            width = dev.width %||% plotdims$W, 
            noRStudioGD = TRUE)
  
  # Set graphical parameters (include mfrow to ensure layout is reverted on exit)
  opar = par(oma = c(0, 0, 3, 0), xpd = NA, mfrow = c(1,1), mar = c(0,0,0,0))
  on.exit(par(opar))
  
  #layout(rbind(1:(1 + length(AM))), widths = widths)
  layout(layoutMat, widths = widths, heights = heights)
  
  # Top margin
  topmar = 0.5 + famnames + !is.null(titles)
  
  # Plot PMs in panel 1 (left)
  if(npm > 0)
    plotPM(PM, nrow = pmDim[1], hatched = hatched, margins = c(2,2,topmar,2), 
           cex = cex, fill = fill, col = col, lwd = lwd, carrier = carrier,...)
  
  # Color style
  miss = dvi$missing
  if(style == 1) 
    fill = list("pink" = miss)
  else if(style == 2) {
    col = list("red" = miss); carrier = miss; lwd = list("1.2" = miss)} 
  
  # Loop through AMs
  nms = names(AM)
  for(i in 1:nam) {
    rw = ceiling(i / amDim[2])
    mar = c(1.5,2,1.5,2)
    if(rw == 1) mar[3] = topmar
    if(rw == amDim[2]) mar[1] = 2
    plot(AM[[i]], hatched = hatched, title = if(famnames) nms[i], margins = mar,
         cex = cex, cex.main = cex+0.2, fill = fill, col = col, lwd = lwd, carrier = carrier, ...)
  }
  
  # X coordinate of each plot region (converted to value in [0,1]).
  ratios = c(0, if(npm > 0) widths[1], sum(widths))/sum(widths)
  
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
  if(npm == 0 && length(titles) == 2)
    titles = titles[2]
  
  midpoints = ratios[1:2] + diff(ratios)/2
  cex.main = args[["cex.main"]] %||% cex + 0.3
  
  font.main = args[["font.main"]] %||% 1
  
  mtext(titles, outer = TRUE, at = midpoints, cex = cex.main, font = font.main)
  
  invisible(c(plotdims, list(layout = layoutMat)))
}

# Default layout of AM families: Up to 5 columns
amArrayDim = function(N, nrow = NA) {
  if(is.na(nrow)) {
    maxcol = if(N <= 9) 4 else if(N <= 15) 5 else if(N<=24) 6 else 7
    nr = (N-1) %/% maxcol + 1
  }
  else 
    nr = min(N, nrow)
  c(nr, ceiling(N/nr))
}

# Default layout of PM singletons
pmArrayDim = function(N, nrow = NA) {
  if(is.na(nrow))
    nr = if(N <= 3) N else ceiling(sqrt(N))
  else 
    nr = min(N, nrow)
  c(nr, ceiling(N/nr))
}


findPlotDims = function(am, amDim = NULL, pmDim = NULL, npm = NULL) {
  amDim = amDim %||% amArrayDim(length(am))
  pmDim = pmDim %||% pmArrayDim(npm)
  
  # Alignment data for each am family
  alignlist = lapply(am, .pedAlignment)
  
  xrange = vapply(alignlist, function(a) max(1, diff(a$xrange)), 1) |> 
    `length<-`(prod(amDim)) |> 
    matrix(nrow = amDim[1], byrow = TRUE)
  
  yrange = vapply(alignlist, function(a) a$maxlev, 1) |> 
    `length<-`(prod(amDim)) |> 
    matrix(nrow = amDim[1], byrow = TRUE)
    
   # AM grid size
  maxWidth = max(rowSums(yrange, na.rm = TRUE))
  maxGen = max(colSums(yrange, na.rm = TRUE))
  
  # Relative widths & heights
  maxXrange = apply(xrange, 2, max, na.rm = TRUE)
  maxYrange = apply(yrange, 1, max, na.rm = TRUE)
  
  widths = sqrt(maxXrange) + 1
  heights = 1 + maxYrange/2
  
  # Add width of PM column if present
  if(pmDim[1] > 0) {
    pmwid = max(sqrt(pmDim[2]), sum(widths)/3)
    widths = c(pmwid, widths)
  }
  
  # Suggested plot dims in inches
  W = sum(widths + 1)
  H = max(3, 1.2 * maxGen) + 0.3
  
  list(amSize = c(maxWidth, maxGen), pmSize = pmDim, 
       widths = widths, heights = heights, W = W, H = H)
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
    densvec = mrg("densvec", 0), 
    fillvec = mrg("fillvec", NA), 
    ltyvec = mrg("ltyvec", 1), 
    lwdvec = mrg("lwdvec", 1),
    carrierTF = mrg("carrierTF", FALSE),
           deceasedTF = mrg("deceasedTF", FALSE))
  
  p = drawPed(align, annot, ...)
  invisible(p)
}



#' Plot DVI solution
#'
#' A version of [plotDVI()] tailor-made to visualise identified individuals, for
#' example as reported by `jointDVI()`.
#'
#' @param dvi A `dviData` object.
#' @param assignment A named character of the format `c(victim = missing, ...)`,
#'   or a data frame produced by [jointDVI()].
#' @param k An integer; the row number when `assignment` is a data frame.
#' @param format A string indicating how identified individuals should be
#'   labelled, using `[M]` and `[S]` as place holders for the missing person and
#'   the matching sample, respectively. (See Examples.)
#' @param ... Further arguments passed on to [plotDVI()].
#'
#' @return NULL.
#'
#' @examples
#'
#' res = jointDVI(example2, verbose = FALSE)
#'
#' plotSolution(example2, res)
#'
#' # With line break in labels
#' plotSolution(example2, res, format = "[M]=\n[S]")
#'
#' # With genotypes for marker 1
#' plotSolution(example2, res, marker = 1)
#'
#' # Non-optimal solutions
#' plotSolution(example2, res, k = 2, pm = FALSE)
#' plotSolution(example2, res, k = 2, cex = 1.3)
#'
#' @export
plotSolution = function(dvi, assignment, k = 1, format = "[S]=[M]", ...) {
  
  # Ensure proper dviData object
  dvi = consolidateDVI(dvi)
  
  a = assignment
  
  if(is.data.frame(a))
    a = unlist(a[k, !names(a) %in% c("loglik", "LR", "posterior"), drop = FALSE])
  
  # Vector of matching pairs
  mtch = a[a != "*"]
  vics = names(mtch)
  
  ### Label format for matching individuals
  fmt = sub("[S]", "%{sample}s", sub("[M]", "%{miss}s", format, fixed = TRUE), fixed = TRUE)
  newlabs = sprintfNamed(fmt, sample = vics, miss = mtch)
  
  # Avoids error if constant format
  if(length(newlabs) != length(mtch))
    newlabs = rep_len(newlabs, length.out = length(mtch))
  
  refs = typedMembers(dvi$am)
  
  if(length(mtch)) {
    am = relabel(dvi$am, old = mtch, new = newlabs)
    am = transferMarkers(dvi$pm[vics], am, idsFrom = vics, idsTo = newlabs, erase = FALSE)
    dvi$am = am
  }
  
  dvi$pm = dvi$pm[a == "*"]
  
  stillmissing = setdiff(dvi$missing, mtch)
  
  plotDVI(dvi, style = 0, fill = list("green" = newlabs), hatched = refs, col = list("red" = stillmissing), 
          lwd = list("1.5" = c(stillmissing)), ...)
}
