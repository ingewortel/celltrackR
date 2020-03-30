#' Track Measures
#'
#' Statistics that can be used to quantify tracks. All of these functions take a single
#' track as input and give a single number as output.
#'
#' @param x a single input track; a matrix whose first column is time and whose
#'  remaining columns are a spatial coordinate.
#' @param from index, or vector of indices, of the first row of the track.
#' @param to index, or vector of indices, of last row of the track.
#' @param xdiff row differences of x.
#' @param degrees logical; should angles be returned in degrees rather than radians?
#'
#' @details
#' Some track measures consider only the first and last position (or steps) of a track,
#' and are most useful in conjunction with \code{\link{aggregate.tracks}}; for instance,
#' \code{squareDisplacement} combined with \code{\link{aggregate.tracks}} gives a mean
#' square displacement plot, and \code{overallAngle} combined with
#' \code{\link{aggregate.tracks}} gives a turning angle plot (see the examples for
#' \code{\link{aggregate.tracks}}). To speed up computation of these measures on
#' subtracks of the same track, the arguments \code{from}, \code{to} and
#' possibly \code{xdiff} are exploited by \code{\link{aggregate.tracks}}.
#'
#' @return
#' \code{trackLength} sums up the distances between subsequent positsion; in other words,
#' it estimates the length of the underlying track by linear interpolation (usually
#' an underestimation). The estimation could be improved in some circumstances by using
#' \code{\link{interpolateTrack}}. The function returns a single, non-negative number.
#'
#' \code{duration} returns the time elapsed between \code{x}'s first and last
#' positions (a single, non-negative number).
#'
#' \code{speed} simply divides \code{\link{trackLength}} by \code{\link{duration}}
#'
#' \code{displacement} returns the Euclidean distance between the track endpoints
#' and \code{squareDisplacement} returns the squared Euclidean distance.
#'
#' \code{displacementVector} returns the vector between the track endpoints. This
#' vector has an element (can be negative) for each (x,y,z) dimension of the coordinates
#' in the track.
#'
#' \code{maxDisplacement} computes the maximal Euclidean distance of any position
#' on the track from the first position.
#'
#' \code{displacementRatio} divides the \code{displacement} by the \code{maxDisplacement};
#' \code{outreachRatio} divides the \code{maxDisplacement} by the \code{trackLength}
#' (Mokhtari et al, 2013). Both measures return
#' values between 0 and 1, where 1 means a perfectly straight track.
#' If the track has \code{trackLength} 0, then \code{NaN} is returned.
#'
#' \code{straightness} divides the \code{displacement} by the \code{trackLength}.
#' This gives a number between 0 and 1, with 1 meaning a perfectly straight track.
#' If the track has \code{trackLength} 0, then \code{NaN} is returned.
#'
#' \code{asphericity} is a different appraoch to measure straightness
#' (Mokhtari et al, 2013): it computes the asphericity of the set of positions on the
#' track _via_ the length of its principal components. Again this gives a number between 0
#' and 1, with higher values indicating straighter tracks.
#' Unlike \code{\link{straightness}}, however, asphericity ignores
#' back-and-forth motion of the object, so something that bounces between two positions
#' will have low \code{straightness} but high \code{asphericity}. We define the
#' asphericity of every track with two or fewer positions to be 1. For one-dimensional
#' tracks with one or more positions, \code{NA} is returned.
#'
#' \code{overallAngle} Computes the angle (in radians) between the first and the last
#' segment of the given track. Angles are measured symmetrically, thus the return values
#' range from 0 to pi; for instance, both a 90 degrees left and right turns yield the
#' value pi/2. This function is useful to generate autocorrelation plots
#' (together with \code{\link{aggregate.tracks}}). Angles can also be returned in degrees,
#' in that case: set \code{degrees = TRUE}.
#'
#' \code{meanTurningAngle} averages the \code{overallAngle} over all
#' adjacent segments of a given track; a low \code{meanTurningAngle} indicates high
#' persistence of orientation, whereas for an uncorrelated random walk we expect
#' 90 degrees. Note that angle measurements will yield \code{NA} values for tracks
#' in which two subsequent positions are identical. By default returns angles in
#' radians; use \code{degrees = TRUE} to return angles in degrees instead.
#'
#' \code{overallDot} computes the dot product between the first and the last
#' segment of the given track. This function is useful to generate autocovariance plots
#' (together with \code{\link{aggregate.tracks}}).
#'
#' \code{overallNormDot} computes the dot product between the unit vectors along
#' the first and the last segment of the given track. This function is useful to
#' generate autocorrelation plots (together with
#' \code{\link{aggregate.tracks}}).
#'
#' \code{hurstExponent} computes the corrected empirical Hurst exponent of the track.
#' This uses the function \code{\link[pracma]{hurstexp}} from the `pracma` package.
#' If the track has less than two positions, NA is returned.
#' \code{fractalDimension} estimates the fractal dimension of a track using the function
#' \code{\link[fractaldim]{fd.estim.boxcount}} from the
#' `fractaldim` package. For self-affine processes in \eqn{n} dimensions,
#' fractal dimension and Hurst exponent
#' are related by the formula \eqn{H=n+1-D}.
#' For non-Brownian motion, however, this relationship
#' need not hold. Intuitively, while the Hurst exponent takes a global approach to the
#' track's properties, fractal dimension is a local approach to the track's properties
#' (Gneiting and Schlather, 2004).
#'
#' @name TrackMeasures
#'
#'
#' @seealso \code{\link{AngleAnalysis}} for methods to compute angles and distances
#'  between pairs of tracks, or of tracks to a reference point, direction, or plane.
#'
#' @examples
#' ## show a turning angle plot with error bars for the T cell data.
#' with( (aggregate(BCells,overallDot,FUN="mean.se",na.rm=TRUE)),{
#'   plot( mean ~ i, xlab="time step",
#'   	ylab="turning angle (rad)", type="l" )
#'   segments( i, lower, y1=upper )
#' } )
#'
#' @references
#' Zeinab Mokhtari, Franziska Mech, Carolin Zitzmann, Mike Hasenberg, Matthias Gunzer
#' and Marc Thilo Figge (2013), Automated Characterization and
#' Parameter--Free Classification of Cell Tracks Based on Local Migration
#' Behavior. \emph{PLoS ONE} \bold{8}(12), e80808. doi:10.1371/journal.pone.0080808
#'
#' Tillmann Gneiting and Martin Schlather (2004), Stochastic Models That Separate Fractal
#' Dimension and the Hurst Effect. \emph{SIAM Review} \bold{46}(2), 269--282.
#' doi:10.1137/S0036144501394387
NULL

#' @rdname TrackMeasures
#' @export
trackLength <- function(x) {
	if (nrow(x) > 2) {
	  # coordinates are in x[,-1,drop=FALSE] (everything but the time col of the track matrix)
	  # apply diff over the columns to get dx,dy,(dz) for every timestep,
	  # compute sqrt( dx^2 + dy^2 + dz^2 ) for every step, then sum up.
		dif <- apply(x[,-1,drop=FALSE], 2, diff)
		return(sum(sqrt(apply(dif^2, 1, sum))))
	} else if (nrow(x) == 2) {
		# this case is necessary because of dimension dropping by 'apply'
		return(sqrt(sum((x[2,-1] - x[1,-1])^2)))

	  # A track of a single coordinate has no length
	} else if (nrow(x) == 1) {
		return(0)

	  # If the matrix is empty, return NA
	} else {
		return(NA)
	}
}

#' @rdname TrackMeasures
#' @export
duration <- function(x) {
  # A track of a single point has a duration of zero
	if( nrow(x) < 2 ){
		return(0)
	}
  # Otherwise, the duration is time difference (col 1)
  # between first and last point
	dur <- x[nrow(x), 1] - x[1,1]
	dur <- unname(dur)
	return(dur)
}

#' @rdname TrackMeasures
#' @export
speed <- function(x) {
  trackLength(x) / duration(x)
}

#' @rdname TrackMeasures
#' @export
displacement <- function( x, from=1, to=nrow(x) ) {
  # tracks of a single coordinate have no displacement
	if( nrow(x) < 2 ){
		return(0)
	}
	sqrt(.rowSums((x[to, -1, drop=FALSE] - x[from, -1, drop=FALSE])^2,
		length(from),ncol(x)-1))
}

#' @rdname TrackMeasures
#' @export
squareDisplacement <- function(x, from=1, to=nrow(x)) {
	if( nrow(x) < 2 ){
		return(0)
	}
	.rowSums((x[to, -1, drop=FALSE] - x[from, -1, drop=FALSE])^2,
		length(from),ncol(x)-1)
}

#' @rdname TrackMeasures
#' @export
displacementVector <- function(x) {
  ret <- x[nrow(x),-1] - x[1,-1]
  rownames(ret) <- NULL
  return(as.vector(ret))
}

#' @rdname TrackMeasures
#' @export
maxDisplacement <- function(x) {
	limits <- c(1,nrow(x))
	sqrt(max(rowSums(sweep(x[seq(limits[1],limits[2]),-1,drop=FALSE],
		2,x[limits[1],-1])^2)))
}

#' @rdname TrackMeasures
#' @export
displacementRatio <- function(x) {
  dmax <- maxDisplacement(x)
  if (dmax > 0) {
    return(displacement(x) / dmax)
  } else {
    return(NaN)
  }
}

#' @rdname TrackMeasures
#' @export
outreachRatio <- function(x) {
  l <- trackLength(x)
  if (l > 0) {
    return(maxDisplacement(x) / l)
  } else {
    return(NaN)
  }
}

#' @rdname TrackMeasures
#' @export
straightness <- function(x) {
  l <- trackLength(x)
  if (l > 0) {
    return(displacement(x) / l)
  } else {
    return(1)
  }
}

#' @rdname TrackMeasures
#' @export
overallAngle <- function (x, from = 1, to = nrow(x), xdiff = diff(x), degrees = FALSE )
{
  # make sure that to and from are the same length
  if( length(from) != length(to) ){
    if( length(from) == 1 ){
      from = rep( from, length(to) )
    } else if ( length(to) == 1 ){
      to = rep( to, length(from) )
    } else {
      stop("overallAngle: from and to must be the same length, or one of them must be a single value." )
    }
  }

  # Start with a zero angle for all angles to compute
  r <- rep(0, length(from))

  # To compute an angle, we need at least two steps, so from must be at least two smaller than to.
  ft <- from < (to - 1)

  # Check if there is anything left to compute
  if (sum(ft) > 0) {
    # Get displacement vectors for the steps to compute angle between.
    # Can be multiple, in a matrix form.
    a <- xdiff[from[ft], -1, drop = FALSE]
    b <- xdiff[to[ft] - 1, -1, drop = FALSE]

    # get angle
    r[ft] <- vecAngle( a, b, degrees = degrees )
  }
  r
}

#' @rdname TrackMeasures
#' @export
meanTurningAngle <- function(x, degrees = FALSE) {
  # angles are only defined for subtracks of at least two steps
  # (three coordinate sets)
	if(nrow(x)<3){
		return(NaN)
	}
  # return the mean angle over all subtracks of 2 steps
	mean(sapply(subtracks(x, 2), overallAngle, degrees = degrees),na.rm=TRUE)
}

#' @rdname TrackMeasures
#' @export
overallDot <- function(x, from=1, to=nrow(x), xdiff=diff(x)) {

  # Start with a NaN angle for all dot products to compute
	r <- rep(NaN, length(from))

	# To compute a dot product, we need at least one step, so from must be at least one smaller than to.
	# all the other ones will return NaN.
	ft <- from<to

	# Check if there is anything left to compute
	if( sum(ft) > 0 ){
	  # Get displacement vectors for the steps to compute angle between.
	  # Can be multiple, in a matrix form.
		a <- xdiff[from[ft],-1,drop=FALSE]
		b <- xdiff[to[ft]-1,-1,drop=FALSE]

		# get dot product
		r[ft] <- .rowSums(a * b, nrow(a), ncol(a))
	}
	r
}

#' @rdname TrackMeasures
#' @export
overallNormDot <- function (x, from = 1, to = nrow(x), xdiff = diff(x))
{
  # Check for tracks of length one
  if(length(xdiff) == 0) { return(NaN) }

  # Drop index column before normalizing
  xdiff <- xdiff[, -1, drop = F]

  # Normalize
  if (requireNamespace("wordspace", quietly = TRUE)) {
    # If wordspace is installed, use faster CPP code
    xdiff_norm <- wordspace::normalize.rows(xdiff)
  } else {
    xdiff_norm <- xdiff / sqrt(rowSums(xdiff^2))
  }

  # In case we're checking for subtracks.length = 1 we know it should be 1
  if(all(from + 1 == to)) {
    return(rep(1, length(from)))
  }

  r <- rep(NaN, length(from))
  ft <- from < to
  if (sum(ft) > 0) {
    a <- xdiff_norm[from[ft], , drop = F]
    b <- xdiff_norm[to[ft] - 1, , drop = F]
    r[ft] <- .rowSums(a * b, nrow(a), ncol(a))
  }
  r
}

#' @rdname TrackMeasures
#' @export
asphericity <- function(x) {
	dim <- ncol(x) - 1
	limits <- c(1,nrow(x))
	if (limits[2]-limits[1]<2) {
		return(1)
	}
	if( dim == 1 ){
		return(NaN)
	}
	eigen.values <- eigen(stats::cov(x[limits[1]:limits[2],-1]))$values
	rav <- mean(eigen.values)
	res <- sum((eigen.values - rav )^2 / dim / (dim - 1) / rav^2)
	return(res)
}

#' @rdname TrackMeasures
#' @export
hurstExponent <- function(x) {
	if( !requireNamespace("pracma",quietly=TRUE) ){
		stop("This function requires the 'pracma' package.")
	}
  if (nrow(x) < 2) {
    return(NA)
  }
  return(pracma::hurstexp(diff(x[,-1,drop=FALSE]), display=FALSE)$Hal)
}

#' @rdname TrackMeasures
#' @export
fractalDimension <- function(x){
  if( !requireNamespace("fractaldim",quietly=TRUE) ){
    stop("This function requires the 'fractaldim' package.")
  }
  return(fractaldim::fd.estim.boxcount(x[,-1])$fd)
}
