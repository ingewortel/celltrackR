#' Plot Tracks in 2D
#'
#' Plots tracks contained in a "tracks" object into a twodimensional space
#' pallelel to the data's axes.
#'
#' @param x the tracks to be plotted.
#' @param dims a vector giving the dimensions of the track data that shall be
#' plotted, e.g. \code{c('x','y')} for the \eqn{x} and \eqn{y} dimension.
#' @param add boolean value indicating whether the tracks are to be added to the
#' current plot.
#' @param col a specification of the color(s) to be used. This can be a vector
#'  of size \code{length(x)}, where each entry specififes the color for the
#'  corresponding track.
#' @param pch.start point symbol with which to label the first position of the track
#'  (see \code{\link[graphics]{points}}).
#' @param pch.end point symbol with which to label the last position of the track
#' @param cex point size for positions on the tracks.
#' @param ... additional parameters (e.g. xlab, ylab).
#' to be passed to \code{\link[graphics]{plot}}
#' (for \code{add=FALSE}) or \code{\link[graphics]{points}} (for \code{add=TRUE}),
#' respectively.
#' @details One dimension of the data (by default \eqn{y}) is plotted against
#' another (by default \eqn{x}). The dimesions can be chosen by means of the
#' parameter \code{dims} and the axes can be labeled accordingly with the aid
#' of \code{xlab} and \code{ylab}. The color can be set through \code{col}.
#' If the tracks should be added to an existing plot, \code{add} is to be set
#' to \code{TRUE}.
#' @seealso \code{\link{plot3d}}
#'
#' @return None
#'
#' @export
plot.tracks <- function(x, dims=c('x','y'), add=F,
                        col=order(names(x)), pch.start=1, pch.end=NULL,
                        cex=.5, ... ) {
  args <- list(...)

  ids <- names(x)
  lcol <- rep_len(col, length(ids))
  pcol <- rep(lcol, sapply(x,nrow))
  dim1 <- unlist( lapply(x,'[', , dims[1]) )
  dim2 <- unlist( lapply(x,'[', , dims[2]) )

  if (add==T) {
    graphics::points( dim1, dim2, col=pcol, cex=cex, ... )
  } else {
    graphics::plot( dim1, dim2, col=pcol, cex=cex, ...)
  }

  if( !is.null(pch.start) ){
    starting.points <- lapply( x, function(p) p[1,] )
    px <- sapply(starting.points,'[[',dims[1])
    py <- sapply(starting.points,'[[',dims[2])
    graphics::points( px, py, col=col, pch=pch.start, cex=2 )
  }

  if( !is.null(pch.end) ){
    end.points <- lapply( x, function(p) p[nrow(p),] )
    px <- sapply(end.points,'[[',dims[1])
    py <- sapply(end.points,'[[',dims[2])
    graphics::points( px, py, col=col, pch=pch.end, cex=2 )
  }

  lseg <- function(f,i) {
    function(d) f( d[,dims[i]], -1 )
  }

  x0 <- .ulapply(x,lseg(utils::head,1))
  y0 <- .ulapply(x,lseg(utils::head,2))

  x1 <- .ulapply(x,lseg(utils::tail,1))
  y1 <- .ulapply(x,lseg(utils::tail,2))

  lcol <- rep(lcol, sapply(x, function(d) {
    nrow(d) - 1
  }))
  graphics::segments(x0, y0, x1, y1, col=lcol)
}


#' Bivariate Scatterplot of Track Measures
#'
#' Plots the values of two measures applied on the given tracks against each
#' other.
#'
#' @param measure.x the measure to be shown on the X axis (see \link{TrackMeasures}).
#' @param measure.y the measure to be shown on the Y axis.
#' @param x the input \code{tracks} object.
#' @param add a logical indicating whether the tracks are to be added to an
#' existing plot via \code{\link[graphics]{points}}.
#' @param xlab label of the x-axis. By default the name of the input function
#' \code{measure.x}.
#' @param ylab label of the y-axis. By default the name of the input function
#' \code{measure.y}.
#' @param ... additional parameters to be passed to \code{\link[graphics]{plot}}
#' (in case \code{add=FALSE}) or \code{\link[graphics]{points}} (\code{add=TRUE}).
# For ellipse documentation:
#' @inheritParams hotellingsTest
#'
#' @details Plots the value of \code{measurey} applied to \code{x} against the
#' value of \code{measurey} applied to \code{y}. This is useful for "FACS-like"
#' motility analysis, where clusters of cell tracks are identified based on their
#' motility parameters (Moreau et al, 2012; Textor et al, 2014).
#'
#' @references
#' Moreau HD, Lemaitre F, Terriac E, Azar G, Piel M, Lennon-Dumenil AM,
#' Bousso P (2012), Dynamic In Situ Cytometry Uncovers
#' T Cell Receptor Signaling during Immunological Synapses and Kinapses In Vivo.
#' \emph{Immunity} \bold{37}(2), 351--363. doi:10.1016/j.immuni.2012.05.014
#'
#' Johannes Textor, Sarah E. Henrickson, Judith N. Mandl, Ulrich H. von Andrian,
#' J\"urgen Westermann, Rob J. de Boer and Joost B. Beltman (2014),
#' Random Migration and Signal Integration Promote Rapid and Robust T Cell Recruitment.
#' \emph{PLoS Computational Biology} \bold{10}(8), e1003752.
#' doi:10.1371/journal.pcbi.1003752
#'
#' @examples
#' ## Compare speed and straightness of 3 example population tracks.
#' ## To make the comparison fair, analyze subtracks of fixed length.
#' plotTrackMeasures( subtracks(TCells,4,0), speed, straightness, ellipse.col="black" )
#' plotTrackMeasures( subtracks(BCells,4,0), speed, straightness,
#'   col=2, ellipse.col=2, pch=2, add=TRUE )
#' plotTrackMeasures( subtracks(Neutrophils,4,0), speed, straightness,
#'   col=3, ellipse.col=3, pch=3, add=TRUE )
#'
#' @return None
#'
#' @export
plotTrackMeasures <- function(x, measure.x, measure.y, add=FALSE,
                              xlab=deparse(substitute(measure.x)),
                              ylab=deparse(substitute(measure.y)),
                              ellipse.col="red", ellipse.border="black",
                              conf.level=0.95, ...) {
  if( !is.tracks(x) ){
    x <- as.tracks( x )
  }
  mx <- sapply(x,measure.x)
  my <- sapply(x,measure.y)
  if(add) {
    graphics::points(mx, my, ...)
  } else {
    graphics::plot(sapply(x,measure.x), sapply(x,measure.y), xlab=xlab, ylab=ylab, ...)
  }
  if( !is.na(ellipse.col) ){
    n <- length(mx)
    p <- 1
    df.1 <- p
    df.2 <- n-p
    d <- cbind( mx, my )
    S <- stats::cov( d )
    t <- sqrt(((n-1)*df.1/(n*df.2))*stats::qf(1-conf.level,df.1,df.2,lower.tail=F))
    graphics::polygon( ellipse::ellipse( S, centre=colMeans(d),
                                         t=t ),
                       border=ellipse.border,
                       col=.setColAlpha(ellipse.col,128) )
  }
}


#' Plot Tracks in 3D
#'
#' Takes an input tracks object and plots them in 3D using the
#' \link[scatterplot3d]{scatterplot3d} function.
#'
#' @param x the tracks which will be plotted in 3d
#' @param ... further arguments to be passed on to
#' \link[scatterplot3d]{scatterplot3d}
#'
#' @examples
#' if( require("scatterplot3d",quietly=TRUE) ){
#'   plot3d( TCells )
#' }
#'
#' @return None.
#' 
#' @export
plot3d <- function(x,...){
  if( !requireNamespace("scatterplot3d",quietly=TRUE) ){
    stop("This function requires the package 'scatterplot3d'; please install manually before continuing: install.packages('scatterplot3d')")
  }
  tracks_df <- as.data.frame.tracks(
    lapply( x, function(t) rbind(t,rep(NA,ncol(t)) ) ) )
  s3d <- scatterplot3d::scatterplot3d(tracks_df[,-c(1,2)],
                                      type="n",xlab="X Position",ylab="Y Position",zlab="Z Position",...)
  colvec <- grDevices::rainbow(length(names(x)))[tracks_df[,1]]
  pts <- s3d$xyz.convert( tracks_df[,-c(1,2)] )
  graphics::segments( utils::head(pts$x,-1), utils::head(pts$y,-1), utils::tail(pts$x,-1),
                      utils::tail(pts$y,-1), col=colvec )
}
