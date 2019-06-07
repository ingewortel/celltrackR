#' Decompose Track(s) into Subtracks
#'
#' Creates a \emph{tracks} object consisting of all subtracks of `x`
#' with `i` segments (i.e., `i`+1 positions).
#'
#' @param x a single track or a \code{tracks} object.
#' @param i subtrack length. A single integer, lists are not supported.
#' @param overlap the number of segments in which each subtrack shall overlap
#' with the previous and next subtrack. The default \code{i - 1} returns all
#' subtracks. Can be negative, which means that space will be left between
#' subtracks.
#'
#' @details The output is always a single \emph{tracks} object, which is
#' convenient for many common analyses. If subtracks are to be considered separately
#' for each track, use the function \code{\link{staggered}} together with
#' \code{lapply}. Subtrack extraction always starts at the first position of the
#' input track.
#'
#' @return A \emph{tracks} object is returned which contains all the subtracks
#' of any track in the input \emph{tracks} object that consist of exactly `i`
#' segments and overlap adjacent subtracks in `overlap` segments.
#'
#' @seealso \code{\link{prefixes}} to extract all subtracks of a given length starting
#' from the first coordinate in each track, \code{\link{subtracksByTime}} to extract
#' all subtracks of a given length  starting at some fixed timepoint,
#' and \code{\link{selectSteps}} to extract single steps starting at a fixed timepoint
#' from a subset of trackids.
#'
#' @export
subtracks <- function(x, i, overlap=i-1 ) {
  if( !is.tracks(x) ){
    x <- wrapTrack( x )
  }
  do.call(c, lapply(x,
                    function(t) .computeSubtracksOfISegments(t, i, overlap )))
}

#' Get Track Prefixes
#'
#' Creates a \code{tracks} object consisting of all prefixes (i.e., subtracks
#' starting with the first position of a track) of `x`
#' with `i` segments (i.e., `i`+1 positions).
#'
#' @param x a single track or a \code{tracks} object.
#' @param i subtrack length. A single integer, lists are not supported.
#'
#' @details This function behaves exactly like \code{\link{subtracks}} except
#' that only subtracks starting from the first position are considered.
#'
#' @return A \emph{tracks} object is returned which contains all the subtracks
#' of any track in the input \emph{tracks} object that consist of exactly `i`
#' segments and start at the first registered coordinate of the given track.
#'
#' @seealso \code{\link{subtracks}} to extract all subtracks of a given length,
#' \code{\link{subtracksByTime}} to extract all subtracks of a given length
#' starting at some fixed timepoint, and \code{\link{selectSteps}} to extract
#' single steps starting at a fixed timepoint from a subset of trackids.
#'
#' @export
prefixes <- function(x,i) {
  if( !is.tracks(x) ){
    x <- wrapTrack( x )
  }
  lapply( x,
          function(t) {
            if( nrow(t) <= i ){ return( NULL ) }
            t[1:(i+1), ,drop=FALSE]
          }
  )
}


#' Extract Subtracks Starting at a Specific Time
#'
#' Obtain all subtracks of i steps (i+1 positions) starting at a given timepoint t.
#'
#' @param X Tracks object to obtain subtracks from.
#' @param t Timepoint at which the subtracks should start. This value is ignored if
#' tlo and thi are specified, see below.
#' @param i Subtrack length (in number of steps). Set this to NULL to obtain subtracks
#' of varying length but within a specified interval [tlo, thi] (see below).
#' @param epsilon Small error allowed when comparing timepoints because of numerical
#' inaccuracies, see details. Timepoints in tracks are returned if they are within
#' [tlo-epsilon, thi+epsilon].
#' @param tlo,thi Interval specifying the timepoints to be returned. By default, these
#' are not used and tracks starting at timepoint t with exactly i steps are returned;
#' see details.
#'
#' @details
#' If i is specified, the given t is retrieved for all tracks in X that contain that
#' timepoint, and any subtracks starting from that time that have exactly i steps are
#' returned. For numerical reasons, timepoints in the data are allowed to deviate
#' a small amount epsilon from t (because otherwise, equal timepoints can seem unequal
#' because of very small deviations).
#'
#' If i is set to NULL, subtracks are returned with all timepoints lying in the interval
#' [ tlo - epsilon, thi + epsilon ]. These subtracks do NOT have to be of equal length.
#'
#'
#' @return A \emph{tracks} object is returned which contains all the subtracks
#' of any track in the input \emph{tracks} object that consist of exactly `i`
#' segments and start at the given timepoint t, OR a
#' \emph{tracks} object with all the timepoints
#' of any track in the input \emph{tracks} object that are between tlo and thi.
#'
#' @seealso \code{\link{subtracks}} to extract all subtracks of a given length,
#' \code{\link{prefixes}} to extract all subtracks of a given length starting
#' from the first coordinate in each track, \code{\link{selectSteps}} to extract
#' single steps starting at a fixed timepoint from a subset of trackids, and
#' \code{\link{timePoints}} to return all timepoints occurring in the dataset.
#'
#' @examples
#' ## Get all the single steps (i=1) starting at the third timepoint in the T cell tracks.
#' subT <- subtracksByTime( TCells, timePoints(TCells)[3], 1 )
#'
#' ## These all have the same number of steps:
#' sapply( subT, nrow )
#'
#' ## Or set i to NULL and return all subtracks within the five first timepoints:
#' subT2 <- subtracksByTime( TCells, NULL, i = NULL,
#'   tlo = timePoints( TCells )[1], thi = timePoints( TCells )[5] )
#'
#' ## These are not all the same length:
#' sapply( subT2, nrow )
#'
#' @export
subtracksByTime <- function( X, t, i = 1, epsilon=1e-4, tlo = t, thi = t )
{
  # Adjust thi and tlo with epsilon
  thi <- thi + epsilon
  tlo <- tlo - epsilon

  # Remove tracks that start at timepoints later than upper bound thi
  # or tracks that end at timepoints before tlo
  X <- X[ sapply(X, function(x) x[1,1] <= thi )]
  if( length(X) != 0 ){
    X <- X[ sapply(X, function(x) x[nrow(x),1] >= tlo )]
  }

  # Filter out parts of the tracks at timepoints starting from tlo
  if( length(X) == 0 ){
    return(X)
  }
  X2 <- as.tracks( lapply( X, function(x) x[ x[,1] >= tlo, , drop=FALSE ] ) )

  if( length(X) == 0 ){
    return(X)
  }

  # If i is specified: get subtracks of length i starting at t by using prefixes
  if( !is.null(i) ){
    sub <- prefixes( X2, i )

  # otherwise: return all subtracks with times between tlo and thi. These do not
  # have to be equal length.
  } else {
    # Filter out parts of the tracks with timepoints until thi
    sub <- lapply( X2, function(x) x[ x[,1] <= thi, , drop=FALSE ] )
  }

  if( length(sub) == 0 ){
    return(sub)
  }

  # Remove any cells for which there are no timepoints
  sub <- sub[ !sapply(sub, is.null) ]

  return(as.tracks(sub))

}


#' Get Single Steps Starting at a Specific Time from a Subset of Tracks
#'
#' Obtain all single steps starting at a given timepoint t from a subset of tracks of interest.
#'
#' @param X Tracks object to obtain subtracks from
#' @param trackids Character vector with the ids of tracks of interest
#' @param t Timepoint at which the subtracks should start
#'
#'
#' @return A \emph{tracks} object is returned which contains all the extracted steps.
#'
#' @seealso \code{\link{subtracks}} to extract all subtracks of a given length,
#' \code{\link{prefixes}} to extract all subtracks of a given length starting
#' from the first coordinate in each track, \code{\link{subtracksByTime}} to
#' extract all subtracks of a given length starting at some fixed timepoint, and
#' \code{\link{timePoints}} to return all timepoints occurring in the dataset.
#'
#' @examples
#' ## Get and plot all steps starting at the third timepoint in tracks 1 and 3 of
#' ## the T cell dataset
#' subT <- selectSteps( TCells, c("1","5"), timePoints(TCells)[3] )
#' plot( subT )
#'
#' @export
selectSteps <- function( X, trackids, t )
{
  X2 <- X[ as.character(trackids) ]
  return( subtracksByTime( X2, t, 1 ) )

}
