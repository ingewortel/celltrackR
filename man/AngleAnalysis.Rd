% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/angle-functions.R
\name{AngleAnalysis}
\alias{AngleAnalysis}
\title{Angle Analysis}
\value{
This page is for documentation only and provides an overview of angle analysis functions
and their use cases. The return values of each of these functions are documented separately;
please follow the link to the documentation page of that specific function.
}
\description{
Analyzing angles to reference directions, points, or planes can be useful to detect
artefacts and/or directionality in tracking datasets (Beltman et al, 2009). All these
functions take a track and a reference (point/direction/plane)
as input and return a distance or angle as output. Angles/distances are by default
computed to the first step in the given track.
}
\details{
\code{\link{angleToPoint}} and \code{\link{distanceToPoint}} return the angle/distance of the track to the
 reference point. The distance returned is between the first coordinate in the track and the
 reference point. The angle is between the overall displacement vector of the track and the
 vector from its first coordinate to the reference point. Angles are by default returned in
 degrees, use \code{degrees=FALSE} to obtain radians. These functions are useful to detect
 directional bias towards a point of interest, which would result in an average angle of less
 than 90 degrees with the reference point (especially for tracks at a small distance to the
 reference point).

\code{\link{angleToPlane}} and \code{\link{distanceToPlane}} return the angle/distance of the track to a
 plane of interest. This plane must be specified by three points lying on it.
 The distance returned is between the first coordinate in the track and the
 reference point. The angle is between the overall displacement vector of the track and the
 plane of interest. These functions are useful to detect tracking artefacts near the borders
 of the imaging volume. Use \code{\link{boundingBox}} to guess where those borders are.
 Angles are by default returned in
 degrees, use \code{degrees=FALSE} to obtain radians.

\code{\link{angleToDir}} returns the angle of a track's overall displacement vector 
 to a direction of interest.
 This function is useful to detect directionality in cases where the direction of the bias is
 known in advance (e.g. when cells are known to move up a chemotactic gradient): in that case,
 the average angle to the reference direction should be less than 90 degrees. Angles are
 by default returned in degrees, use \code{degrees=FALSE} to obtain radians.

\code{\link{angleSteps}} and \code{\link{distanceSteps}} return the angle/distance between a pair
 of steps in the data that occur at the same timepoint. Angles are in degrees by default,
 use \code{degrees=FALSE} to obtain radians. Use \code{\link{stepPairs}} to extract all pairs of
 steps that occur at the same timepoint, and use \code{\link{analyzeStepPairs}} to do this and then
 also obtain the angles and distances for each of these pairs.

\code{\link{angleCells}} and \code{\link{distanceCells}} return the angle/distance between a pair
 of tracks in the data. The computed angles are between the overall displacement vectors of the
 tracks, the distance is the shortest distance between them at any timepoint they share.
 Angles are in degrees by default, use \code{degrees=FALSE} to obtain radians.
 Use \code{\link{cellPairs}} to extract all pairs of
 cells in the data, and use \code{\link{analyzeCellPairs}} to do this and then
 also obtain the angles and distances for each of these pairs.
}
\examples{
## Plotting the angle versus the distance to a reference point can be informative to
## detect biased movement towards that point. We should be suspicious especially
## when small angles are more frequent at lower distances.
steps <- subtracks( sample( Neutrophils, 50 ), 1 )
bb <- boundingBox( Neutrophils )
angles <- sapply( steps, angleToPoint, p = bb["max",-1] )
distances <- sapply( steps, distanceToPoint, p = bb["max",-1] )
scatter.smooth( distances, angles )
abline( h = 90, col = "red" )

## Get a distribution of Neutrophil step angles with the reference direction
## in positive y direction. The histogram is enriched for low angles, suggesting
## directed movement:
hist( sapply( steps, angleToDir, dvec=c(1,-1) ) )

## Plotting the angle versus the distance to a reference plane can be informative to
## detect tracking artefacts near the border of the imaging volume.
## We should be suspicious especially when small angles are more frequent at low distances
## to the border planes; as is the case in the z-dimension for the raw data:
load( system.file("extdata", "TCellsRaw.rda", package="celltrackR" ) )
steps <- subtracks( sample( TCellsRaw, 50 ), 1 )
minz <- boundingBox( TCellsRaw )["min","z"]
## Compute angles and distances to the lower plane in z-dimension
angles <- sapply( steps, angleToPlane, p1 = c(0,0,minz), p2 = c(1,0,minz), p3 = c(0,1,minz) )
distances <- sapply( steps, distanceToPlane, p1 = c(0,0,minz), p2 = c(1,0,minz), p3 = c(0,1,minz) )
scatter.smooth( distances, angles )
abline( h = 32.7, col = "red" )

## Plot distance versus angle for all cell pairs (here in only a sample to speed things up)
pairs <- analyzeCellPairs( sample( TCells, 50 ) )
scatter.smooth( pairs$dist, pairs$angle )
abline( h = 90, col = "red" )

## Plot distance versus angle for all step pairs, filtering for those that
## displace at least 2 microns
pairs <- analyzeStepPairs( sample( TCells, 50 ), filter.steps = function(t) displacement(t) > 2 )
scatter.smooth( pairs$dist, pairs$angle )
abline( h = 90, col = "red" )


}
\references{
Joost B. Beltman, Athanasius F.M. Maree and Rob. J. de Boer (2009),
Analysing immune cell migration. \emph{Nature Reviews Immunology} \bold{9},
789--798. doi:10.1038/nri2638
}
\seealso{
\code{\link{TrackMeasures}} for other measures that can be used to quantify tracks.

See the vignettes on Quality Control and Track Analysis for more detailed examples of
angle analyses.
\code{browseVignettes( package = "celltrackR" )}
}
