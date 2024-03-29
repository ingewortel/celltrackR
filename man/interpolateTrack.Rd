% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{interpolateTrack}
\alias{interpolateTrack}
\title{Interpolate Track Positions}
\usage{
interpolateTrack(x, t, how = "linear")
}
\arguments{
\item{x}{the input track (a matrix or data frame).}

\item{t}{the times at which to approximate track positions. These must lie
within the interval spanned by the track timepoints.}

\item{how}{specifies how to perform the interpolation. Possible values are
\code{"linear"} (which uses \code{\link[stats:approxfun]{approx}} with default values) and
\code{"spline"} (which uses \code{\link[stats:splinefun]{spline}} with default values).}
}
\value{
The interpolated track (a matrix or data frame).
}
\description{
Approximates the track positions at given time points using linear interpolation
(via the \code{\link[stats:approxfun]{approx}} function).
}
\examples{
## Compare interpolated and non-interpolated versions of a track
bb <- boundingBox( TCells[2] )
plot( TCells[2] )
t2i <- interpolateTrack(TCells[[2]], seq(bb[1,"t"],bb[2,"t"],length.out=100),"spline")
plot( tracks( t2i ), add=TRUE, col=2 )

}
