#' Filter Tracks
#'
#' Extracts subtracks based on a given function.
#'
#' @param f a function that accepts a single track as its first argument and returns a
#' logical value (or a value that can be coerced to a locical).
#' @param x a tracks object.
#' @param ... further arguments to be passed on to \code{f}.
#' @return A \code{tracks} object containing only those tracks from \code{x} for which
#' \code{f} evaluates to \code{TRUE}.
#' @examples
#' ## Remove short tracks from the T cells data
#' plot( filterTracks( function(t) nrow(t)>10, TCells ) )
#' @export
filterTracks <- function(f,x,...){
	structure(x[as.logical(sapply(x,function(x) f(x,...) ))], class="tracks")
}

#' Sort Track Positions by Time
#'
#' Sorts the positions in each track in a \emph{tracks} object by time.
#'
#' @param x the \emph{tracks} object whose tracks are to be sorted by time.
#'
#' @param decreasing logical.  Should the sort be increasing or decreasing?
#'  Provided only for consistency with the generic sort method. The positions in
#'  each track should be sorted in increasing time order.
#' @param ... further arguments to be passed on to \code{order}.
#'
#' @details Sorts the positions of each track (represented as a data frame) in the
#' \emph{tracks} object by time (given in the column t).
#'
#' @return A \emph{tracks object} that contains the tracks from the input object
#' sorted by time is returned.
#' @export
sort.tracks <- function(x, decreasing=FALSE, ...) {
	as.tracks(lapply(x, function(d) d[order(d[,"t"],decreasing=decreasing),,drop=FALSE] ))
}


#' Read Tracks from Text File
#'
#' Reads cell tracks from a CSV or other text file. Data are expected to be organized as
#' follows.
#' One column contains a track identifier, which can be numeric or a string, and
#' determines which points belong to the same track.
#' Another column is expected to contain a time index or a time period (e.g. number of
#' seconds elapsed since the beginning of the track, or since the beginning of the
#' experiment). Input of dates is not (yet) supported, as absolute time information is
#' frequently not available.
#' Further columns contain the spatial coordinates. If there are three or less spatial
#' coordinates, their names will by "x", "y", and "z"
#' (depending on whether the tracks are 1D, 2D or 3D). If there are four or more spatial
#' coordinates, their names will be "x1", "x2", and so on.
#' The names or indices of these columns in the CSV files are given using the
#' corresponding parameters (see below). Names and indices can be mixed, e.g. you can
#' specify \code{id.column="Parent"} and \code{pos.columns=1:3}
#'
#' @param file the name of the file which the data are to be read from, a
#' readable text-mode connection or a complete URL
#' (see \code{\link[utils]{read.table}}).
#'
#' @param id.column index or name of the column that contains the track ID.
#'
#' @param time.column index or name of the column that contains elapsed time.
#'
#' @param pos.columns vector containing indices or names of the columns that contain
#' the spatial coordinates. If this vector has two entries and the second entry is NA,
#' e.g. \code{c('x',NA)} or \code{c(5,NA)} then all columns from the indicated column
#' to the last column are used. This is useful when reading files where the exact number
#' of spatial dimensions is not known beforehand.
#'
#' @param scale.t a value by which to multiply each time point. Useful for changing units,
#' or for specifying the time between positions if this is not contained in the file
#' itself.
#'
#' @param scale.pos a value, or a vector of values, by which to multiply each spatial
#' position. Useful for changing units.
#'
#' @param header a logical value indicating whether the file contains the
#' names of the variables as its first line. See \code{\link[utils]{read.table}}.
#'
#' @param sep a character specifying how the colums of the data are separated.
#' The default value \code{""} means columns are separated by tabs or other spaces.
#' See \code{\link[utils]{read.table}}.
#'
#' @param track.sep.blankline logical. If set to \code{TRUE}, then tracks are expected
#' to be separated by one or more blank lines in the input file instead of being
#' designated by a track ID column. In this case, numerical track IDs are automatically
#' generated.
#'
#' @param ... further arguments to be passed to \code{read.csv}, for instance
#' \code{sep="\t"} can be useful for tab-separated files.
#'
#' @return An object of class \emph{tracks} is returned, which is a list of
#' matrices, each containing the positions of one track. The matrices
#' have a column \eqn{t}, followed by one column for each of the input track's
#' coordinates.
#'
#' @details The input file's first four fields are interpreted as \eqn{id},
#' \eqn{pos}, \eqn{t} and \eqn{x}, respectively, and, if available, the fifth
#' as \eqn{y} and the sixth as \eqn{z}. The returned object has the class
#' \emph{tracks}, which is a list of data frames representing the single
#' tracks and having columns \eqn{t} and \eqn{x}, plus \eqn{y} and \eqn{z}, if
#' necessary. The tracks' ids are retained in their position in the list, while
#' the field \eqn{pos} will be unmaintained.
#'
#' @export
read.tracks.csv <- function(file, id.column=1, time.column=2,
	pos.columns=c(3,4,5), scale.t=1, scale.pos=1,
	header=TRUE, sep="", track.sep.blankline=FALSE, ...) {
	if( track.sep.blankline ){
		data.raw <- utils::read.table(file, header=header, sep=sep,
			blank.lines.skip = FALSE, ...)
		is.blank <- apply( data.raw, 1, function(x) all( is.na(x) ) )
		ids <- cumsum( is.blank )
		data.raw <- cbind( data.raw, ids )
		data.raw <- data.raw[!is.blank,]
		id.column <- ncol( data.raw )
	} else {
		data.raw <- utils::read.table(file, header=header, sep=sep, ...)
	}
	as.tracks.data.frame(data.raw, id.column, time.column, pos.columns, scale.t, scale.pos)
}

#' Split Track into Multiple Tracks
#'
#' @param x the input track (a data frame or a matrix).
#' @param id a string used to identify the resulting tracks; for instance,
#' if \code{id="X"}, then the resulting tracks are named X_1, X_2 and so forth.
#' Otherwise, they are simply labelled with integer numbers.
#' @param positions a vector of positive integers, given in ascending order.
#' @param min.length nonnegative integer. Resulting tracks that have fewer positions than
#'  the value of this parameter are dropped.
#'
#' @return An object of class \emph{tracks} with the resulting splitted tracks.
#' 
#' @export
splitTrack <- function( x, positions, id=NULL, min.length=2 ){
	freqs <- diff(c(0,positions,nrow(x)))
	segs <- rep( seq_len(length(positions)+1), freqs )
	segs[rep(freqs,freqs) < min.length] <- NA
	r <- structure(
		split.data.frame( as.matrix(x), segs ),
		class="tracks")
	if(!is.null(id) && ( length(r) > 0 ) ){
		if( length(r) > 1 ){
			names( r ) <- paste0( id, "_", seq_along(r) )
		} else if( length(r) == 1 ){
			names( r ) <- id
		}
	}
	return(r)
}




#' Staggered Version of a Function
#'
#' Returns the "staggered" version of a track measure. That is, instead of
#' computing the measure on the whole track, the measure is averaged over
#' all subtracks (of any length) of the track.
#'
#' @param measure a track measure (see \link{TrackMeasures}).
#' @param ... further parameters passed on to \code{\link{applyStaggered}}.
#'
#' @details This is a wrapper mainly designed to provide a convenient interface
#' for track-based staggered computations with \code{lapply}, see example.
#'
#' @return Returns a function that computes
#' the given measure in a staggered fashion on that track.
#'
#' @references
#' Zeinab Mokhtari, Franziska Mech, Carolin Zitzmann, Mike Hasenberg, Matthias Gunzer
#' and Marc Thilo Figge (2013), Automated Characterization and
#' Parameter--Free Classification of Cell Tracks Based on Local Migration
#' Behavior. \emph{PLoS ONE} \bold{8}(12), e80808. doi:10.1371/journal.pone.0080808
#'
#' @examples
#' hist( sapply( TCells, staggered( displacement ) ) )
#'
#' @export
staggered <- function(measure, ...){
	return( function(x) applyStaggered(x, measure, ...) )
}


#' Compute a Measure on a Track in a Staggered Fashion
#'
#' Computes a measure on all subtracks of a track and return them either
#' as a matrix or return their mean.
#'
#' @param x the track for which the measure is to be computed.
#' @param measure the measure that is to be computed.
#' @param matrix a logical indicating whether the whole matrix of values for
#' the measure for each of the input track's subtracks is to be returned.
#' Otherwise only the mean is returned.
#' @param min.segments the number of segments that each regarded subtrack
#' should at least consist of. Typically, this value would be set to the
#' minimum number of segments that a (sub)track must have in order for the
#' measure to be decently computed. For example, at least two segments are needed
#' to compute the \code{\link{overallAngle}}.
#'
#' @details The measure is computed for each of the input track's subtracks of
#' length at least \code{min.segments}, and the resulting values are either
#' returned in a matrix (if \code{matrix} is set), or their mean is returned.
#' The computed matrix is symmetric since the direction along which a
#' track is traversed is assumed not to matter. The values at
#' \code{[i, i + j]}, where j is a nonnegative integer with
#' \eqn{j < }\code{min.segments}, (with the default value \code{min.segments=1}
#' this is exactly the main diagonal) are set to \code{NA}.
#'
#' @return If \code{matrix} is set, a matrix with the values of the measure for
#' all the input track's subtracks is returned. The value of this matrix at
#' position \code{[i, j]} corresponds to the subtrack that starts with the input track's
#' \eqn{j}th point and ends at its \eqn{i}th. Thus, with increasing column number,
#' the regarded subtrack's starting point is advanced on the original track,
#' and the values for increasingly long subtracks starting from this point can
#' be found columnwise below the main diagonal, respectively.
#' If `matrix` is not set, the mean over the values of the measure for all
#' subtracks of at least `min.segments` segments is retruned.
#'
#' @examples
#' ## Compute the staggered matrix for overallAngle applied to all long enough
#' ## subtracks of the first T cell track
#' applyStaggered(TCells[[1]], overallAngle, matrix=TRUE, min.segments = 2)
#' @export
applyStaggered <- function(x, measure, matrix=FALSE, min.segments=1) {
	if (matrix) {
		mat <- matrix(nrow=nrow(x), ncol=nrow(x), 0)
		diag(mat) <- NA
	} else {
		stag.meas <- c()
	}
	if (!matrix) {
		segMeans <- .computeSegmentwiseMeans(x, measure, min.segments)
		n.before <- length(segMeans)
		segMeans <- segMeans[!is.nan(segMeans)]
		n <- length(segMeans)
		if (n.before!=n) {
		  warning("NaNs have been removed.")
		}
		weights <- seq(n, 1)
		weights <- weights[!is.nan(segMeans)]
		ret <- stats::weighted.mean(segMeans, weights)
		return(ret)
	} else {
		for (i in (min.segments):(nrow(x) - 2)) {
			subtracks <- subtracks(x, i)
			val <- sapply(subtracks, measure)
			diag(mat[-1:-i,]) <- val
		}
		mat[nrow(x), 1] <- sapply(subtracks(x, (nrow(x)-1)), measure)
		mat <- mat + t(mat)
		return(mat)
	}
}


#' Normalize a Measure to Track Duration
#'
#' Returns a measure that divides the input measure by the duration of its
#' input track.
#'
#' @param x a track measure (see \link{TrackMeasures}).
#'
#' @return A function that computes the input measure for a given track
#' and returns the result divided by the track's duration.
#'
#' @examples
#' ## normalizeToDuration(displacement) can be used as an indicator
#' ## for the motion's efficiency
#' sapply(TCells, normalizeToDuration(displacement))
#' @export
normalizeToDuration <- function(x) {
  function(track) {
    x(track) / duration(track)
  }
}










#' Compute Summary Statistics of Subtracks
#'
#' Computes a given measure on subtracks of a given track set, applies a summary
#' statistic for each subtrack length, and returns the results in a convenient form.
#' This important workhorse function facilitates many common motility analyses
#' such as mean square displacement, turning angle, and autocorrelation plots.
#'
#' @param x the tracks object whose subtracks are to be considered.
#' If a single track is given, it will be coerced to a tracks object
#' using \code{\link{wrapTrack}} (but note that this requires an explicit call
#' \code{aggregate.tracks}).
#' @param measure the measure that is to be computed on the subtracks.
#' @param by a string that indicates how grouping is performed. Currently, two
#' kinds of grouping are supported:
#' \itemize{
#'  \item "subtracks"  Apply \code{measure} to all subtracks according to
#' the parameters \code{subtrack.length} and \code{max.overlap}.
#'  \item "prefixes"   Apply \code{measure} to all prefixes (i.e., subtracks starting
#'  from a track's initial position) according to the parameter \code{subtrack.length}.
#' }
#'
#' @param FUN a summary statistic to be computed on the measures of subtracks with the
#' same length. Can be a function or a string.
#' If a string is given, it is first matched to the following builtin values:
#' \itemize{
#'  \item "mean.sd"  Outputs the mean and \eqn{mean - sd} as lower and
#'  \eqn{mean + sd} as upper bound 
#'  \item "mean.se" Outputs the mean and \eqn{mean - se} as lower and
#'  \eqn{mean + se} as upper bound 
#'  \item "mean.ci.95" Outputs the mean and upper and lower bound of a
#'  parametric 95 percent confidence interval.
#'  \item "mean.ci.99" Outputs the mean and upper and lower bound of a
#'  parametric 95 percent confidence intervall.
#'  \item "iqr" Outputs the interquartile range, that is, the median, and the
#'  25-percent-quartile as a lower and and the 75-percent-quartile as an
#'  upper bound
#' }
#' If the string is not equal to any of these, it is passed on to
#' \code{\link{match.fun}}.
#'
#' @param subtrack.length an integer or a vector of integers defining which subtrack
#' lengths are considered. In particular, \code{subtrack.length=1}
#' corresponds to a "step-based analysis" (Beltman et al, 2009).
#'
#' @param max.overlap an integer controlling what to do with overlapping subtracks.
#' A maximum overlap of \code{max(subtrack.length)} will imply
#' that all subtracks are considered. For a maximum overlap of 0, only non-overlapping
#' subtracks are considered. A negative overlap can be used to ensure that only subtracks
#' a certain distance apart are considered. In general, for non-Brownian motion there will
#' be correlations between subsequent steps, such that a negative overlap may be necessary
#' to get a proper error estimate.
#' @param na.rm logical. If \code{TRUE}, then \code{NA}'s and \code{NaN}'s
#' are removed prior to computing the summary statistic.
#' @param filter.subtracks a function that can be supplied to exclude certain subtracks
#' from an analysis. For instance, one may wish to compute angles only between steps of
#' a certain minimum length (see Examples).
#' @param count.subtracks logical. If \code{TRUE}, the returned dataframe contains an
#' extra column \code{ntracks} showing the number of subtracks of each length. This is 
#' useful to keep track of since the returned \code{value} estimates for high 
#' \code{i} are often based on very few subtracks.
#' @param ... further arguments passed to or used by methods.
#'
#' @details For every number of segments \eqn{i} in the set defined by
#' \code{subtrack.length}, all subtracks of any track in the input
#' \code{tracks} object that consist of exactly \eqn{i} segments are
#' considered. The input \code{measure} is applied to the subtracks individually,
#' and the \code{statistic} is applied to the resulting values.
#'
#' @return A data frame with one row for every \eqn{i}
#' specified by \code{subtrack.length}. The first column contains the values
#' of \eqn{i} and the remaining columns contain the values of the summary statistic
#' of the measure values of tracks having exactly \eqn{i} segments.
#' @examples
#' ## A mean square displacement plot with error bars.
#' dat <- aggregate(TCells, squareDisplacement, FUN="mean.se")
#' with( dat ,{
#'   plot( mean ~ i, xlab="time step",
#'   	ylab="mean square displacement", type="l" )
#'   segments( i, lower, y1=upper )
#' } )
#' 
#' ## Note that the values at high i are often based on very few subtracks:
#' msd <- aggregate( TCells, squareDisplacement, count.subtracks = TRUE )
#' tail( msd )
#'
#' ## Compute a turning angle plot for the B cell data, taking only steps of at least
#' ## 1 micrometer length into account
#' check <- function(x) all( sapply( list(head(x,2),tail(x,2)), trackLength ) >= 1.0 )
#' plot( aggregate( BCells, overallAngle, subtrack.length=1:10,
#'   filter.subtracks=check )[,2], type='l' )
#'
#' ## Compare 3 different variants of a mean displacement plot
#' # 1. average over all subtracks
#' plot( aggregate( TCells, displacement ), type='l' )
#' # 2. average over all non-overlapping subtracks
#' lines( aggregate( TCells, displacement, max.overlap=0 ), col=2 )
#' # 3. average over all subtracks starting at 1st position
#' lines( aggregate( TCells, displacement, by="prefixes" ), col=3 )
#'
#' @references
#' Joost B. Beltman, Athanasius F.M. Maree and Rob. J. de Boer (2009),
#' Analysing immune cell migration. \emph{Nature Reviews Immunology} \bold{9},
#' 789--798. doi:10.1038/nri2638
#'
#' @importFrom stats aggregate
#'
#' @aliases aggregate
#' @export
aggregate.tracks <- function( x, measure, by="subtracks", FUN=mean,
    subtrack.length=seq(1, (maxTrackLength(x)-1)),
    max.overlap=max(subtrack.length), na.rm=FALSE,
    filter.subtracks=NULL, count.subtracks=FALSE, ... ){
    if( !is.tracks( x ) ){
    	if( is.data.frame(x) || is.matrix(x) ){
    		x <- wrapTrack( x )
    	} else {
    		stop("Cannot coerce argument 'tracks' to tracks object" )
    	}
    }
    if (max(subtrack.length) > (maxTrackLength(x)-1)) {
		warning("No track is long enough!")
  	}
  	the.statistic <- NULL
  	if (is.character(FUN)) {
		if (FUN == "mean.ci.95") {
		  the.statistic <- function(x) {
		  	mx <- mean(x,na.rm=na.rm)
			ci <- tryCatch( stats::t.test(x)$conf.int, error=function(e) rep(NA,2) )
			return(c(lower=ci[1], mean=mx, upper=ci[2]))
		  }
		} else if (FUN == "mean.ci.99") {
		  the.statistic <- function(x) {
			mx <- mean(x,na.rm=na.rm)
			ci <- tryCatch( stats::t.test(x, conf.level=.99)$conf.int,
				error=function(e) rep(NA,2) )
			return(c(lower=ci[1], mean=mean(x), upper=ci[2]))
		  }
		} else if (FUN == "mean.se") {
		  the.statistic <- function(x) {
		  	mx <- mean(x,na.rm=na.rm)
		  	sem <- stats::sd(x,na.rm=na.rm)/sqrt(sum(!is.na(x)))
		  		# note that this also works is na.omit is FALSE
			return(c(mean=mx, lower = mx - sem, upper = mx + sem))
		  }
		} else if (FUN == "mean.sd") {
		  the.statistic <- function(x) {
		  	mx <- mean(x,na.rm=na.rm)
		  	sem <- stats::sd(x,na.rm=na.rm)
			return(c(mean=mx, lower = mx - sem, upper = mx + sem))
		  }
		} else if (FUN == "iqr") {
		  the.statistic <- function(x) {
			lx <- stats::quantile(x,probs=c(.25,.5,.75),na.rm=na.rm)
			names(lx) <- c("lower","median","upper")
			return(lx)
		  }
		} else {
			the.statistic.raw <- match.fun( FUN )
		}
	} else {
		if ( !is.function( FUN) ) {
		   stop("FUN must be a function or a string")
		} else {
			the.statistic.raw <- FUN
		}
	}
	if( is.null( the.statistic ) ){
		if( na.rm ){
			the.statistic <- function(x){
				the.statistic.raw(x[!is.na(x) & !is.nan(x)])
			}
		} else {
			the.statistic <- the.statistic.raw
		}
	}
	if( is.function( filter.subtracks ) ){
		isValidSubtrack <- function(t) !is.null(t) && filter.subtracks(t)
	} else {
		isValidSubtrack <- Negate(is.null)
	}
	if( by == "subtracks" ){
		# optimized version: avoids making copies of subtracks (relevant especially
		# for measures like displacement, which actually only consider the track
		# endpoints). Also pre-computes the diff of the track if needed.
		if( ( length( intersect(c("x","limits"),names(formals(measure))) ) == 2 )
			&& !is.function( filter.subtracks ) ){
			if( "xdiff" %in% names(formals(measure)) ){
				ft <- function(t,i) apply( .subtrackIndices(i,min(max.overlap,i-1),nrow(t)),
					1, measure, x=t, xdiff=diff(t) )
			} else {
				ft <- function(t,i) apply( .subtrackIndices(i,min(max.overlap,i-1),nrow(t)),
					1, measure, x=t )
			}
			measure.values <- list()
			k <- 1
			xi <- x
			for( i in sort(subtrack.length) ){
				xi <- xi[sapply(xi,nrow)>i]
				measure.values[[k]] <- .ulapply( xi, ft, i=i )
				k <- k+1
			}
		} else if( ( length( intersect(c("x","from","to"),names(formals(measure))) ) == 3 )
			&& !is.function( filter.subtracks ) ){
			measure.values <- list()
			k <- 1
			xi <- x
			for( i in sort(subtrack.length) ){
				xi <- xi[sapply(xi,nrow)>i]
				xi1 <- do.call( rbind, xi )
				si <- .subtrackIndices(i,min(max.overlap,i-1),sapply(xi,nrow))
				measure.values[[k]] <- measure(xi1,si[,1],si[,2])
				k <- k+1
			}
		} else {
			# unoptimized version: makes copies of all subtracks.
			the.subtracks <- list()
			measure.values <- list()
			k <- 1
			for (i in subtrack.length) {
				the.subtracks[[k]] <- subtracks(x, i, min(max.overlap,i-1) )
				the.subtracks[[k]] <- Filter( isValidSubtrack, the.subtracks[[k]] )
				measure.values[[k]] <- sapply( the.subtracks[[k]], measure )
				k <- k+1
			}
		}
	} else if (by == "prefixes" ) {
		measure.values <- lapply( subtrack.length,
				function(i) sapply( Filter( isValidSubtrack, prefixes( x, i )  ),
					measure )
			)
	} else {
		stop( "grouping method ", by, " not supported" )
	}
	# this is necessary because of automatic insertion of NULL elements when e.g.
	# measure.values[[5]] is assigned to
	measure.values <- Filter( Negate(is.null), measure.values )
	value <- sapply(measure.values, the.statistic)	
	ret <- rbind(i=subtrack.length, value)
	if(count.subtracks){
		num <- sapply( measure.values, length )
		ret <- rbind( ret, ntracks = num )
	}
	return(data.frame(t(ret)))
}

#' Length of Longest Track
#'
#' Determines the maximum number of positions over the tracks in \code{x}.
#'
#' @param x the \code{tracks} object the tracks in which are to be considered.
#'
#' @return The maximum number of rows of a track in \code{x}
#'
#' @export
maxTrackLength <- function(x) {
  return(max(sapply(x, nrow)))
}


#' Bounding Box of a Tracks Object
#'
#' Computes the minimum and maximum coordinates per dimension (including time)
#' for all positions in a given list of tracks.
#'
#' @param x the input \code{tracks} object.
#'
#' @return Returns a matrix with two rows and \eqn{d+1} columns, where \eqn{d} is
#' the number of spatial dimensions of the tracks. The first row contains the minimum
#' and the second row the maximum value of any track in the dimension given by
#' the column.
#'
#' @examples
#' ## Use bounding box to set up plot window
#' bb <- boundingBox(c(TCells,BCells,Neutrophils))
#' plot( Neutrophils, xlim=bb[,"x"], ylim=bb[,"y"], col=1 )
#' plot( BCells, col=2, add=TRUE )
#' plot( TCells, col=3, add=TRUE )
#'
#' @export
boundingBox <- function(x) {
	if( !is.tracks(x) ){
		x <- as.tracks(x)
	}
	bounding.boxes <- lapply(x, .boundingBoxTrack)
	empty <- mat.or.vec(nrow(bounding.boxes[[1]]), ncol(bounding.boxes[[1]]))
	colnames(empty) <- colnames(bounding.boxes[[1]])
	rownames(empty) <- rownames(bounding.boxes[[1]])
	empty[1,] <- Inf
	empty[2,] <- -Inf
	return(Reduce(.minMaxOfTwoMatrices, bounding.boxes, empty))
}

#' Test Unbiasedness of Motion
#'
#' Test the null hypothesis that a given set of tracks originates from an uncorrelated
#' and unbiased type of motion (e.g., a random walk without drift). This is done by
#' testing whether the mean step vector is equal to the null vector.
#'
#' @param tracks the tracks whose biasedness is to be determined.
#' @param dim vector with the names of the track's
#' dimensions that are to be considered. By default c("x", "y").
#' @param plot logical indicating whether the scatter of the step's directions,
#' origin of ordinates (green circle) and the mean of the data points (green
#' cross) are to be plotted. (In one dimension also the bounds of the
#' condfidence interval are given.) Plot works only in one or two dimensions.
#' @param add whether to add the plot to the current plot (\code{TRUE}) or create a
#" new one (\code{FALSE}).
#' @param step.spacing How many positions are to be left out between
#' the steps that are considered for the test. For persistent motion, subsequent
#' steps will be correlated, which leads to too low p-values because Hotelling's
#' test assumes that the input data is independent. To avoid this, the resulting
#' p-value should either be corrected for this dependence (e.g. by adjusting
#' the degrees of freedom accordingly), or `step.spacing` should be set to a value
#' high enough to ensure that the considered steps are approximately independent.
#' @param ellipse.col color with which to draw the confidence ellipse of the mean (for
#'  1D, this corresponds to the confidence interval of the mean).
#'  Use \code{NA} to omit the confidence ellipse.
#' @param ellipse.border color of the confidence ellipse border. Use \code{NA} to omit
#'  the border.
#' @param conf.level the desired confidence level for the confidence ellipse.
#' @param ... further arguments passed on to \code{plot}.
#'
#' @return A list with class \code{htest}.
#'
#' @details Computes the displacement vectors of all segments in the tracks
#' given in \code{tracks}, and performs Hotelling's T-square Test on that vector.
#'
#' @examples
#' ## Test H_0: T-cells migrate by uncorrelated random walk on x and y coordinates,
#' ## and report the p-value.
#' hotellingsTest( TCells )$p.value
#'
#' @references
#' Johannes Textor, Antonio Peixoto, Sarah E. Henrickson, Mathieu
#'  Sinn, Ulrich H. von Andrian and Juergen Westermann (2011),
#'	Defining the Quantitative Limits of Intravital Two-Photon Lymphocyte Tracking.
#' \emph{PNAS} \bold{108}(30):12401--12406. doi:10.1073/pnas.1102288108
#'
#' @export
hotellingsTest <- function(tracks, dim=c("x", "y"),
	step.spacing=0, plot=FALSE, add=FALSE, ellipse.col="blue",
	ellipse.border="black", conf.level=0.95, ... ) {
	if( plot && !(length(dim) %in% c(1,2) ) ){
		stop("Plotting is only supported for 1D and 2D")
	}
	Tx <- projectDimensions(tracks, dim)
	sx <- sapply( subtracks(Tx, 1, overlap=-step.spacing), displacementVector )
	if( length(dim) > 1 ){
		sx <- t(sx)
	} else {
		sx <- matrix(sx,ncol=1)
	}
	n <- dim(sx)[1]
	p <- dim(sx)[2]
	df.1 <- p
	df.2 <- n-p
	S <- stats::cov(sx)
	mX <- colMeans(sx)
	A <- solve(S)
	T2 <- n * t(mX) %*% A %*% mX
	p.value <- 1-stats::pf(T2*df.2/df.1/(n-1), df.1, df.2)
	if (plot && (length(dim)==2)) {
		if (plot) {
			if( add ){
				graphics::points(sx,...)
			} else {
				graphics::plot(sx,...)
			}
			graphics::abline(h=0)
			graphics::abline(v=0)
			if( !is.na(ellipse.col) ){
				t <- sqrt(((n-1)*df.1/(n*df.2))*stats::qf(1-conf.level,df.1,df.2,lower.tail=F))
				graphics::polygon(ellipse::ellipse(S,centre=mX,t=t),
					col=.setColAlpha(ellipse.col,128),border=ellipse.border)
				graphics::points(mX[1],mX[2],col=ellipse.col,pch=20)
			}
		}
	} else if (plot && (length(dim)==1)) {
		sd.sx <- stats::sd(sx)
		graphics::plot(0, 0, col=3, xlab = dim[1], ylim=c(-1,1), ...)
		graphics::stripchart(sx, add=TRUE, at=0, pch=1)
		graphics::abline(v=0)
		t <- sqrt(((n-1)*df.1/(n*df.2))*stats::qf(1-conf.level,df.1,df.2,lower.tail=F))
		graphics::polygon(ellipse::ellipse(rbind(c(S,0),c(0,.1)),centre=c(mX,0),t=t),
			col=.setColAlpha(ellipse.col,128),border=NA)
		graphics::points(mX[1],0,col=ellipse.col,pch=20)
	}
	structure(list(statistic=c(T2=T2),
		p.value=p.value,
		method="Hotelling's one sample T2-test",
		parameter=c(df1=df.1,df2=df.2),
		alternative="two.sided",
		null.value=c(location=paste0("c(", paste0(rep(0,ncol(sx)), collapse = ","), ")"))
		),
		class="htest")
}

#' Select Tracks by Measure Values
#'
#' Given a tracks object, extract a subset based on upper and lower bounds of a certain
#' measure. For instance, extract all tracks with a certain minimum length.
#'
#' @param x the input tracks.
#' @param measure measure on which the selection is based (see \link{TrackMeasures}).
#' @param lower specifies the lower bound (inclusive) of the allowable measure.
#' @param upper specifies the upper bound (inclusive) of the allowable measure.
#'
#' @return A \code{\link{tracks}} object with the selected subset of tracks.
#' @examples
#' ## Slower half of T cells
#' slow.tcells <- selectTracks( TCells, speed, -Inf, median( sapply(TCells,speed) ) )
#'
#' @export
selectTracks <- function(x,measure,lower,upper){
	if(is.tracks(x)){
		mdat <- sapply(x,measure)
		return(x[mdat<=upper & mdat>=lower])
	}else{
		obj <-deparse(substitute(tracks))
		warning(sprintf("Obj:%s is not a tracks object",obj))
	}
}


#' Compute Time Step of Tracks
#'
#' Applies a summary statistics on the time intervals between pairs of consecutive
#' positions in a track dataset.
#'
#' @param x the input tracks.
#' @param FUN the summary statistic to be applied.
#' @param na.rm logical, indicates whether to remove missing values before applying FUN.
#'
#' @details Most track quantification depends on the assumption that track positions are
#' recorded at constant time intervals. If this is not the case, then most of the statistics
#' in this package (except for some very simple ones like \code{\link{duration}}) will not work.
#' In reality, at least small fluctuations of the time steps can be expected. This
#' function provides a means for quality control with respect to the tracking time.
#'
#' @return Summary statistic of the time intervals between two consecutive positions in a track
#' dataset.
#'
#' @examples
#' ## Show tracking time fluctuations for the T cell data
#' d <- timeStep( TCells )
#' plot( sapply( subtracks( TCells, 1 ) , duration ) - d, ylim=c(-d,d) )
#'
#' @importFrom stats median
#' @export
timeStep <- function( x, FUN=median, na.rm=FALSE ){
	if( is.character(FUN) ){
		FUN=match.fun( FUN )
	}
	x <- .ulapply( x, function(x) diff(x[,1]) )
	if( na.rm ){
		x <- x[!is.na(x)]
	}
	return( FUN( x ) )
}

#' Interpolate Track Positions
#'
#' Approximates the track positions at given time points using linear interpolation
#' (via the \code{\link[stats:approxfun]{approx}} function).
#'
#' @param x the input track (a matrix or data frame).
#' @param t the times at which to approximate track positions. These must lie
#' within the interval spanned by the track timepoints.
#' @param how specifies how to perform the interpolation. Possible values are
#' \code{"linear"} (which uses \code{\link[stats:approxfun]{approx}} with default values) and
#' \code{"spline"} (which uses \code{\link[stats:splinefun]{spline}} with default values).
#'
#' @examples
#' ## Compare interpolated and non-interpolated versions of a track
#' bb <- boundingBox( TCells[2] )
#' plot( TCells[2] )
#' t2i <- interpolateTrack(TCells[[2]], seq(bb[1,"t"],bb[2,"t"],length.out=100),"spline")
#' plot( tracks( t2i ), add=TRUE, col=2 )
#'
#' @return The interpolated track (a matrix or data frame).
#'
#' @export
interpolateTrack <- function( x, t, how="linear" ){
	f=stats::approx
	if( how=="spline" ){
		f=stats::spline
	}
	pos.interpolated <- apply(x[,-1],2,function(p) f(x[,1], p, xout=t)$y)
	r <- cbind( t, pos.interpolated )
	colnames(r) <- colnames(x)
	return(r)
}


#' Subsample Track by Constant Factor
#'
#' Make tracks more coarse-grained by keeping only every \emph{k}th position.
#'
#' @param x an input track or tracks object.
#' @param k a positive integer. Every \eqn{k}th position of each
#' input track is kept.
#' @seealso \code{interpolateTrack}, which can be used for more flexible track
#' coarse-graining.
#' @return A \code{\link{tracks}} object with the new, more coarse-grained tracks.
#' @examples
#' ## Compare original and subsampled versions of the T cell tracks
#' plot( TCells, col=1 )
#' plot( subsample( TCells, 3 ), col=2, add=TRUE, pch.start=NULL )
#'
#' @export
subsample <- function( x, k=2 ){
	if( is.tracks(x) ){
		as.tracks( lapply( x, function(t) subsample(t,k) ) )
	} else {
		x[seq(1,nrow(x),by=k),]
	}
}


#' Find All Unique Time Points in a Track Dataset
#'
#' Return a vector of all the timepoints t found in any of the matrices in the
#' tracks object.
#'
#' @param X a tracks object.
#'
#' @return A numeric vector of unique timepoints.
#'
#'
#' @examples
#' ## Get all timepoints in the T cell dataset
#' tp <- timePoints( TCells )
#'
#' @export
timePoints <- function( X )
{
  timepoints.per.track <- lapply( X, function(x) unique( x[,1] ) )
  return( sort( unique( unlist( timepoints.per.track ) ) ) )
}
