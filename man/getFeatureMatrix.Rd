% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster-functions.R
\name{getFeatureMatrix}
\alias{getFeatureMatrix}
\title{Obtaining A Feature Matrix}
\usage{
getFeatureMatrix(tracks, measures, dist = FALSE, ...)
}
\arguments{
\item{tracks}{the tracks that are to be analyzed.}

\item{measures}{a function, or a vector of functions (see \link{TrackMeasures}).
Each function is expected to return a single number given a single track.}

\item{dist}{should a distance matrix rather than a feature matrix be returned?}

\item{...}{further arguments passed on to "dist"}
}
\value{
A matrix with a row for each track and a column for each measure.
}
\description{
Applies a given vector of track measures directly on a set of tracks, returning
output in a matrix with a column for each measure and a row for each track. Can
also return a distance matrix, which some clustering methods require.
}
\examples{
## Get speed, meanTurningAngle, and straightness for T cell tracks
fm <- getFeatureMatrix( TCells, c(speed,meanTurningAngle,straightness))
str(fm)

}
\seealso{
\code{\link{clusterTracks}} for a quick method to compute the feature
matrix and a clustering, and \code{\link{trackFeatureMap}} to perform
dimensionality reduction methods on a set of track features.
}
