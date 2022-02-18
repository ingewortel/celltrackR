#' Angle Analysis
#'
#' Analyzing angles to reference directions, points, or planes can be useful to detect
#' artefacts and/or directionality in tracking datasets (Beltman et al, 2009). All these
#' functions take a track and a reference (point/direction/plane)
#' as input and return a distance or angle as output. Angles/distances are by default
#' computed to the first step in the given track.
#'
#'
#' @details
#'
#' \code{\link{angleToPoint}} and \code{\link{distanceToPoint}} return the angle/distance of the track to the
#'  reference point. The distance returned is between the first coordinate in the track and the
#'  reference point. The angle is between the overall displacement vector of the track and the
#'  vector from its first coordinate to the reference point. Angles are by default returned in
#'  degrees, use \code{degrees=FALSE} to obtain radians. These functions are useful to detect
#'  directional bias towards a point of interest, which would result in an average angle of less
#'  than 90 degrees with the reference point (especially for tracks at a small distance to the
#'  reference point).
#'
#' \code{\link{angleToPlane}} and \code{\link{distanceToPlane}} return the angle/distance of the track to a
#'  plane of interest. This plane must be specified by three points lying on it.
#'  The distance returned is between the first coordinate in the track and the
#'  reference point. The angle is between the overall displacement vector of the track and the
#'  plane of interest. These functions are useful to detect tracking artefacts near the borders
#'  of the imaging volume. Use \code{\link{boundingBox}} to guess where those borders are.
#'  Angles are by default returned in
#'  degrees, use \code{degrees=FALSE} to obtain radians.
#'
#' \code{\link{angleToDir}} returns the angle of a track's overall displacement vector 
#'  to a direction of interest.
#'  This function is useful to detect directionality in cases where the direction of the bias is
#'  known in advance (e.g. when cells are known to move up a chemotactic gradient): in that case,
#'  the average angle to the reference direction should be less than 90 degrees. Angles are
#'  by default returned in degrees, use \code{degrees=FALSE} to obtain radians.
#'
#' \code{\link{angleSteps}} and \code{\link{distanceSteps}} return the angle/distance between a pair
#'  of steps in the data that occur at the same timepoint. Angles are in degrees by default,
#'  use \code{degrees=FALSE} to obtain radians. Use \code{\link{stepPairs}} to extract all pairs of
#'  steps that occur at the same timepoint, and use \code{\link{analyzeStepPairs}} to do this and then
#'  also obtain the angles and distances for each of these pairs.
#'
#' \code{\link{angleCells}} and \code{\link{distanceCells}} return the angle/distance between a pair
#'  of tracks in the data. The computed angles are between the overall displacement vectors of the
#'  tracks, the distance is the shortest distance between them at any timepoint they share.
#'  Angles are in degrees by default, use \code{degrees=FALSE} to obtain radians.
#'  Use \code{\link{cellPairs}} to extract all pairs of
#'  cells in the data, and use \code{\link{analyzeCellPairs}} to do this and then
#'  also obtain the angles and distances for each of these pairs.
#'
#' @name AngleAnalysis
#'
#' @seealso \code{\link{TrackMeasures}} for other measures that can be used to quantify tracks.
#'
#' See the vignettes on Quality Control and Track Analysis for more detailed examples of
#' angle analyses.
#' \code{browseVignettes( package = "celltrackR" )}
#'
#' @return This page is for documentation only and provides an overview of angle analysis functions
#' and their use cases. The return values of each of these functions are documented separately;
#' please follow the link to the documentation page of that specific function.
#'
#' @examples
#' ## Plotting the angle versus the distance to a reference point can be informative to
#' ## detect biased movement towards that point. We should be suspicious especially
#' ## when small angles are more frequent at lower distances.
#' steps <- subtracks( sample( Neutrophils, 50 ), 1 )
#' bb <- boundingBox( Neutrophils )
#' angles <- sapply( steps, angleToPoint, p = bb["max",-1] )
#' distances <- sapply( steps, distanceToPoint, p = bb["max",-1] )
#' scatter.smooth( distances, angles )
#' abline( h = 90, col = "red" )
#'
#' ## Get a distribution of Neutrophil step angles with the reference direction
#' ## in positive y direction. The histogram is enriched for low angles, suggesting
#' ## directed movement:
#' hist( sapply( steps, angleToDir, dvec=c(1,-1) ) )
#'
#' ## Plotting the angle versus the distance to a reference plane can be informative to
#' ## detect tracking artefacts near the border of the imaging volume.
#' ## We should be suspicious especially when small angles are more frequent at low distances
#' ## to the border planes; as is the case in the z-dimension for the raw data:
#' load( system.file("extdata", "TCellsRaw.rda", package="celltrackR" ) )
#' steps <- subtracks( sample( TCellsRaw, 50 ), 1 )
#' minz <- boundingBox( TCellsRaw )["min","z"]
#' ## Compute angles and distances to the lower plane in z-dimension
#' angles <- sapply( steps, angleToPlane, p1 = c(0,0,minz), p2 = c(1,0,minz), p3 = c(0,1,minz) )
#' distances <- sapply( steps, distanceToPlane, p1 = c(0,0,minz), p2 = c(1,0,minz), p3 = c(0,1,minz) )
#' scatter.smooth( distances, angles )
#' abline( h = 32.7, col = "red" )
#'
#' ## Plot distance versus angle for all cell pairs (here in only a sample to speed things up)
#' pairs <- analyzeCellPairs( sample( TCells, 50 ) )
#' scatter.smooth( pairs$dist, pairs$angle )
#' abline( h = 90, col = "red" )
#'
#' ## Plot distance versus angle for all step pairs, filtering for those that
#' ## displace at least 2 microns
#' pairs <- analyzeStepPairs( sample( TCells, 50 ), filter.steps = function(t) displacement(t) > 2 )
#' scatter.smooth( pairs$dist, pairs$angle )
#' abline( h = 90, col = "red" )
#'
#'
#' @references
#' Joost B. Beltman, Athanasius F.M. Maree and Rob. J. de Boer (2009),
#' Analysing immune cell migration. \emph{Nature Reviews Immunology} \bold{9},
#' 789--798. doi:10.1038/nri2638
NULL


#' Angle Between Two Vectors
#'
#' Compute the angle between two vectors a and b, which can be numeric vectors
#' or matrices in which each row represents a  numeric vector.
#' In the last case, one angle is returned for each row. By default, angles
#' are returned in degrees -- set \code{degrees = TRUE} to return radians.
#'
#' @param a the first vector or set of vectors. Must be a numeric vector or a
#'  matrix where each row represents a numeric vector.
#' @param b the second vector or set of vectors, for which angles with the
#'  vector (set) a must be computed. Must have the same dimensions as a.
#' @param degrees logical: if \code{TRUE} (default), return angles in degrees instead of radians.
#'
#' @return A single angle (if a and b are single vectors) or a numeric vector of
#'  angles (if a and b are matrices; in that case, the output vector contains one
#'  angle for each row in matrices a and b).
#'
#' @examples
#' ## The angle between the vectors [0,1] and [1,0] is 90 degrees:
#' vecAngle( c(0,1), c(1,0) )
#' ## The same holds for 3D angles:
#' vecAngle( c(0,1,0), c(1,0,0) )
#' @export
vecAngle <- function( a, b, degrees = TRUE )
{
  # Check if inputs have the right dimensions
  if( !( is.numeric(a) && is.numeric(b) && (is.matrix(a) == is.matrix(b)) ) ){
    stop( "vecAngle: a and b must both be numeric vectors of the same length or matrices of the same size")
  }
  if( any( pracma::size(a) != pracma::size(b) ) ){
    stop( "vecAngle: cannot compute angle between vectors of unequal dimensions.")
  }

  # Convert vector to matrix if necessary
  if( !is.matrix(a) ){
    a <- matrix( a, ncol = length(a) )
    b <- matrix( b, ncol = length(b) )
  }

  # Normalize a and b by their length
  a <- a/sqrt(.rowSums(a^2, nrow(a), ncol(a)))
  b <- b/sqrt(.rowSums(b^2, nrow(b), ncol(b)))

  # Compute dot product
  rs <- .rowSums(a * b, nrow(a), ncol(a))
  rs[rs > 1] <- 1
  rs[rs < -1] <- -1

  # angle is the acos of the dot product of the normalized vectors
  ang <- acos(rs)

  if( degrees ){
    return( pracma::rad2deg(ang) )
  } else{
    return(ang)
  }
}

#' Angle with a Reference Point
#'
#' Compute the angle between a track's overall displacement vector and the vector from
#' it's first coordinate to a reference point. Useful to
#' detect directed movement towards a point (see examples).
#'
#' @param x a single input track; a matrix whose first column is time and whose
#'  remaining columns are a spatial coordinate.
#' @param p numeric vector of coordinates of the reference point p to compute angles/distances to.
#' @param from index, or vector of indices, of the first row of the track. If
#' \code{from} is a vector, angles are returned for all steps starting at
#' the indices in \code{from}.
#' @param degrees logical; should angles be returned in degrees rather than radians? (default = TRUE).
#'
#'
#' @return A single angle.
#'
#' @details The average angle of steps to a reference point should be 90 degrees if there is
#'  no bias towards movement in the direction of the reference point. If there is such a bias,
#'  there should be an enrichment of smaller angles. The expected distribution without bias
#'  is a uniform distribution in 2D or a sine distribution in 3D (Beltman et al, 2009).
#'
#' @seealso \code{\link{distanceToPoint}} to compute the distance to the reference point, and
#'  \code{\link{AngleAnalysis}} for other methods to compute angles and distances.
#'
#' @references
#' Joost B. Beltman, Athanasius F.M. Maree and Rob. J. de Boer (2009),
#' Analysing immune cell migration. \emph{Nature Reviews Immunology} \bold{9},
#' 789--798. doi:10.1038/nri2638
#'
#' @examples
#' ## Get a distribution of step angles with a reference point
#' ## Use bb to get the corner with highest x,y (,z) value
#' ## The histogram is enriched for low angles, suggesting directed movement:
#' steps <- subtracks( Neutrophils, 1 )
#' bb <- boundingBox( Neutrophils )
#' hist( sapply( steps, angleToPoint, p = bb["max",-1] ) )
#'
#' ## The same does not hold for movement of T cells towards the point (0,0)
#' steps <- subtracks( TCells, 1 )
#' hist( sapply( steps, angleToPoint, p = c(0,0) ) )
#'
#' ## Plotting the angle versus the distance to the reference point can also be informative,
#' ## especially when small angles are more frequent at lower distances.
#' angles <- sapply( steps, angleToPoint, p = bb["max",-1] )
#' distances <- sapply( steps, distanceToPoint, p = bb["max",-1] )
#' scatter.smooth( distances, angles )
#' abline( h = 90, col = "red" )
#' @export
angleToPoint <- function (x, p = c(1,1,1), from = 1, degrees = TRUE )
{

  # Check if the given point has the correct dimensions
  if( length(p) != ncol(x) - 1 ){
    stop("In angleToPoint: Reference point coordinates must have the same number of dimensions as
         coordinates in the tracking data.")
  }

  # If point is equal to the first point in the track, the angle is undefined.
  if( all( p == x[1,-1] ) ){
    warning("In angleToPoint: Reference point is equal to the track starting point. Angle
            is undefined, returning NA.")
    return(NA)
  }

  # Initialize a vector of NAs for every index in "from".
  # This will be replaced if possible.
  r <- rep(NA, length(from))

  # To compute an angle we need at least one step (two timepoints).
  ft <- from < nrow(x)
  if( sum(ft) > 0 ){
    # x,y,z components of the step vectors a to compute angle with refpoint with
    starts <- x[from[ft], -1, drop = FALSE]
    # displacements from these 'from' points to the track endpoint
    a <- t( apply( starts, 1, function(t) x[nrow(x),-1] - t ) )

    # compute vectors of step starting points with the refpoint
    b <- x[from[ft],-1,drop=FALSE]
    b <- t( apply( b, 1, function(x) p - x ) )

    # compute angle
    r[ft] <- vecAngle( a, b, degrees = degrees )
  }
  r
}


#' Angle with a Reference Direction
#'
#' Compute the angle between a track's overall displacement and a reference direction.
#' Useful to detect biased movement when the directional bias is known (see examples).
#'
#' @param x a single input track; a matrix whose first column is time and whose
#'  remaining columns are a spatial coordinate.
#' @param dvec numeric vector specifying a reference direction to compute angles to.
#' @param from index, or vector of indices, of the first row of the track. If
#' \code{from} is a vector, angles are returned for all steps starting at
#' the indices in \code{from}.
#' @param degrees logical; should angles be returned in degrees rather than radians? (default = TRUE).
#'
#' @return A single angle.
#'
#' @details The average angle of steps to a reference direction should be 90 degrees if there is
#'  no bias towards movement in the direction of the reference point. If there is such a bias,
#'  there should be an enrichment of smaller angles. The expected distribution without bias
#'  is a uniform distribution in 2D or a sine distribution in 3D (Beltman et al, 2009).
#'
#' @seealso \code{\link{AngleAnalysis}} for other methods to compute angles and distances.
#'
#' @references
#' Joost B. Beltman, Athanasius F.M. Maree and Rob. J. de Boer (2009),
#' Analysing immune cell migration. \emph{Nature Reviews Immunology} \bold{9},
#' 789--798. doi:10.1038/nri2638
#'
#' @examples
#' ## Get a distribution of Neutrophil step angles with the reference direction in positive
#' ## y direction. The histogram is enriched for low angles, suggesting directed movement:
#' steps <- subtracks( Neutrophils, 1 )
#' hist( sapply( steps, angleToDir, dvec=c(1,-1) ) )
#' @export
angleToDir <- function (x, dvec = c(1,1,1), from = 1, degrees=TRUE )
{

  # Check if the given direction has the correct dimensions
  if( length(dvec) != ncol(x) - 1 ){
    stop("In angleToDir: Direction vector must have the same number of dimensions as
         coordinates in the tracking data.")
  }

  # Initialize a vector of NAs for every index in "from".
  # This will be replaced if possible.
  r <- rep(NA, length(from))

  # To compute an angle we need at least one step (two timepoints).
  ft <- from < nrow(x)
  if( sum(ft) > 0 ){
    # x,y,z components of the step vectors a to compute angle with dvec with
    starts <- x[from[ft], -1, drop = FALSE]
    # displacements from these 'from' points to the track endpoint
    a <- t( apply( starts, 1, function(t) x[nrow(x),-1] - t ) )
    b <- matrix( rep( dvec, nrow(a)),
                 ncol = length(dvec), byrow = TRUE )

    # Normalize step vectors and the ref direction vector by their length
    #a <- a/sqrt(.rowSums(a^2, nrow(a),ncol(a)))
    #b <- dvec/sqrt( sum( dvec^2 ) )

    # Dot product of normalized a and b; this equals cos alpha
    #rs <- .rowSums(a * b, nrow(a), ncol(a))
    #rs[rs > 1] <- 1
    #rs[rs < -1] <- -1

    # angle is acos.
    r[ft] <- vecAngle( a, b, degrees = degrees )
  }
  r
  }


#' Angle with a Reference Plane
#'
#' Compute the angle between a track's overall displacement and a reference plane.
#' Useful to detect directed movement and/or tracking artefacts.
#'
#' @param x a single input track; a matrix whose first column is time and whose
#'  remaining columns are a spatial coordinate.
#' @param p1,p2,p3 numeric vectors of coordinates of three points specifying a reference plane to
#'  compute distances to.
#' @param from index, or vector of indices, of the first row of the track. If
#' \code{from} is a vector, angles are returned for all steps starting at
#' the indices in \code{from}.
#' @param degrees logical; should angles be returned in degrees rather than radians? (default = TRUE).
#'
#' @return A single angle.
#'
#' @details The average angle of steps to a reference plane should be roughly 32.7 degrees.
#'  Lower angles to the border planes of an imaging volume can be indicative of tracking
#'  artefacts, and systematic deviations from 32.7 can indicate a directional bias
#'  (Beltman et al, 2009).
#'
#' @seealso \code{\link{distanceToPlane}} to compute the distance to the reference plane, and
#'  \code{\link{AngleAnalysis}} for other methods to compute angles and distances.
#'
#' @references
#' Joost B. Beltman, Athanasius F.M. Maree and Rob. J. de Boer (2009),
#' Analysing immune cell migration. \emph{Nature Reviews Immunology} \bold{9},
#' 789--798. doi:10.1038/nri2638
#'
#' @examples
#' ## Plotting the angle versus the distance to a reference plane can be informative to
#' ## detect tracking artefacts near the border of the imaging volume.
#' ## We should be suspicious especially when small angles are more frequent at low distances
#' ## to the border planes.
#' load( system.file("extdata", "TCellsRaw.rda", package="celltrackR" ) )
#' steps <- subtracks( TCellsRaw, 1 )
#' minz <- boundingBox( TCellsRaw )["min","z"]
#' ## Compute angles and distances to the lower plane in z-dimension
#' angles <- sapply( steps, angleToPlane, p1 = c(0,0,minz), p2 = c(1,0,minz), p3 = c(0,1,minz) )
#' distances <- sapply( steps, distanceToPlane, p1 = c(0,0,minz), p2 = c(1,0,minz), p3 = c(0,1,minz) )
#' scatter.smooth( distances, angles )
#' abline( h = 32.7, col = "red" )
#' @export
angleToPlane <- function (x, p1 = c(0,0,0), p2 = c(0,1,0), p3 = c(1,0,0),
                          from = 1, degrees =TRUE )
{

  # Check if the given points have the correct dimensions
  plengths <- c( length(p1), length(p2), length(p3) )
  if( length(unique(plengths)) > 1 ){
    stop("In angleToPlane: Points p1,p2,p3 specifying the plane must have the same number of coordinates.")
  }
  if( unique(plengths) != 3 || ncol(x)-1 != 3 ){
    stop("In angleToPlane: Method is only supported for three-dimensional data.")
  }
  pmatrix <- matrix( c(p1,p2,p3), nrow = 3, byrow = TRUE )
  if( any( duplicated( pmatrix ) ) ){
    stop("In angleToPlane: Points p1, p2, and p3 must be three unique points!")
  }
  if( length(p1) != ( ncol(x) - 1 ) ){
    stop("In angleToPlane: Plane points must have the same number of coordinates as
         coordinates in the tracking data.")
  }

  # Check if the given points actually span a plane and are not on the same line.
  vec1 <- p2 - p1
  vec2 <- p3 - p1
  if( vecAngle( vec1, vec2 ) == 0 || vecAngle( vec1, vec2 ) == 180 ){
    stop("In angleToPlane: Points p1, p2, and p3 are on the same line and do not fully specify a plane.")
  }

  # Initialize a vector of NAs for every index in "from".
  # This will be replaced if possible.
  r <- rep(NA, length(from))

  # To compute an angle we need at least one step (two timepoints).
  ft <- from < nrow(x)
  if( sum(ft) > 0 ){
    # x,y,z components of the step vectors a to compute angle with plane with
    starts <- x[from[ft], -1, drop = FALSE]
    # displacements from these 'from' points to the track endpoint
    a <- t( apply( starts, 1, function(t) x[nrow(x),-1] - t ) )

    # Compute the normal vector of the plane. Use the three points to get
    # two vectors in the plane, then use their cross product to get the normal.
    planematrix <- matrix( c(p1,p2,p3), ncol = length(p1), byrow=TRUE )
    planevecs <- diff(planematrix)
    pnorm <- pracma::crossn( planevecs )
    pnorm <- matrix( rep( pnorm, nrow(a)),
                     ncol = length(pnorm), byrow = TRUE )

    # # Normalize step vectors and the plane normal vector by their length
    # #a <- a/sqrt(.rowSums(a^2, nrow(a),ncol(a)))
    # b <- pnorm/sqrt( sum( pnorm^2 ) )
    #
    # # Dot product of normalized a and b; this equals cos alpha
    # rs <- .rowSums(a * b, nrow(a), ncol(a))
    # rs[rs > 1] <- 1
    # rs[rs < -1] <- -1

    # angle to the plane normal vector
    rtmp <- vecAngle( a, pnorm, degrees = degrees )

    # angle to the plane itself is 90-r (or pi/2 - r )
    deg90 <- ifelse( degrees, 90, pi/2 )
    rnew <- abs( deg90 - rtmp )
    r[ft] <- rnew
  }
  return(r)
  }


#' Distance to a Reference Plane
#'
#' Compute the (shortest) distance between the starting point of a track and a reference plane.
#' Useful to detect directed movement and/or tracking artefacts.
#'
#' @param x a single input track; a matrix whose first column is time and whose
#'  remaining columns are a spatial coordinate.
#' @param p1,p2,p3 numeric vectors of coordinates of three points specifying a reference plane to
#'  compute distances to.
#' @param from index, or vector of indices, of the first row of the track. If
#' \code{from} is a vector, distances are returned for all steps starting at
#' the indices in \code{from}.
#'
#' @return A single distance.
#'
#' @seealso \code{\link{angleToPlane}} to compute the angle to the plane, and
#'  \code{\link{AngleAnalysis}} for other methods to compute angles and distances.
#'
#' @examples
#' ## Plotting the angle versus the distance to a reference plane can be informative to
#' ## detect tracking artefacts near the border of the imaging volume.
#' ## We should be suspicious especially when small angles are more frequent at low distances
#' ## to the border planes.
#' load( system.file("extdata", "TCellsRaw.rda", package="celltrackR" ) )
#' steps <- subtracks( TCellsRaw, 1 )
#' minz <- boundingBox( TCellsRaw )["min","z"]
#' ## Compute angles and distances to the lower plane in z-dimension
#' angles <- sapply( steps, angleToPlane, p1 = c(0,0,minz), p2 = c(1,0,minz), p3 = c(0,1,minz) )
#' distances <- sapply( steps, distanceToPlane, p1 = c(0,0,minz), p2 = c(1,0,minz), p3 = c(0,1,minz) )
#' scatter.smooth( distances, angles )
#' abline( h = 32.7, col = "red" )
#' @export
distanceToPlane <- function (x, p1 = c(0,0,0), p2 = c(0,1,0), p3 = c(1,0,0), from = 1  )
{

  # Check if the given points have the correct dimensions
  plengths <- c( length(p1), length(p2), length(p3) )
  if( length(unique(plengths)) != 1){
    stop("In distanceToPlane: Points p1,p2,p3 specifying the plane must have the same number of coordinates.")
  }
  if( unique(plengths) != 3 || ncol(x)-1 != 3 ){
    stop("In distanceToPlane: Method is only supported for three-dimensional data.")
  }
  pmatrix <- matrix( c(p1,p2,p3), nrow = 3, byrow = TRUE )
  if( any( duplicated( pmatrix ) ) ){
    stop("In distanceToPlane: Points p1, p2, and p3 must be three unique points!")
  }
  if( length(p1) != ( ncol(x) - 1 ) ){
    stop("In distanceToPlane: Plane points must have the same number of coordinates as
         coordinates in the tracking data. Currently, the method only supports 3D data.")
  }

  # Check if the given points actually span a plane and are not on the same line.
  vec1 <- p2 - p1
  vec2 <- p3 - p1
  if( vecAngle( vec1, vec2 ) == 0 || vecAngle( vec1, vec2 ) == 180 ){
    stop("In distanceToPlane: Points p1, p2, and p3 are on the same line and do not fully specify a plane.")
  }

  # Initialize a vector of NAs for every index in "from".
  # This will be replaced if possible.
  r <- rep(NA, length(from))

  # To compute a step distance we need at least one step (two timepoints).
  ft <- from < nrow(x)
  if( sum(ft) > 0 ){
    # x,y,z components of the step startpoints a to compute distance to plane for
    a <- x[from[ft], -1, drop = FALSE]

    # Compute the normal vector of the plane. Use the three points to get
    # two vectors in the plane, then use their cross product to get the normal.
    planematrix <- matrix( c(p1,p2,p3), ncol = length(p1), byrow=TRUE )
    planevecs <- diff(planematrix)
    pnorm <- pracma::crossn( planevecs )
    pnorm.len <- sqrt( sum( pnorm^2 ) )


    # Compute distance:
    D <- -sum( pnorm * p1 )
    distances <- unname( apply( a, 1, function(x)
      abs(sum(x*pnorm)+D)/pnorm.len ) )

    r[ft] <- distances
  }
  r
  }

#' Distance to a Reference Point
#'
#' Compute the distance between the starting point of a track and a reference point.
#' Useful to
#' detect directed movement towards a point (see examples).
#'
#' @param x a single input track; a matrix whose first column is time and whose
#'  remaining columns are a spatial coordinate.
#' @param p numeric vector of coordinates of the reference point p to compute distances to.
#' @param from index, or vector of indices, of the first row of the track. If
#' \code{from} is a vector, distances are returned for all steps starting at
#' the indices in \code{from}.
#'
#' @return A single distance.
#'
#' @seealso \code{\link{angleToPoint}} to compute the angle to the reference point, and
#'  \code{\link{AngleAnalysis}} for other methods to compute angles and distances.
#'
#' @examples
#' ## Plotting the angle versus the distance to a reference point can be informative to
#' ## detect biased movement towards that point. We should be suspicious especially
#' ## when small angles are more frequent at lower distances.
#' steps <- subtracks( Neutrophils, 1 )
#' bb <- boundingBox( Neutrophils )
#' angles <- sapply( steps, angleToPoint, p = bb["max",-1] )
#' distances <- sapply( steps, distanceToPoint, p = bb["max",-1] )
#' scatter.smooth( distances, angles )
#' abline( h = 90, col = "red" )
#' @export
distanceToPoint <- function (x, p = c(0,0,0), from = 1 )
{
  # Check if the given point has the correct dimensions
  if( length(p) != ncol(x) - 1 ){
    stop("In distanceToPoint: Reference point coordinates must have the same number of dimensions as coordinates in the tracking data.")
  }

  # Initialize a vector of NAs for every index in "from".
  # This will be replaced if possible.
  r <- rep(NA, length(from))

  # To compute step distance to point we need at least one step (two timepoints).
  ft <- from < nrow(x)
  if( sum(ft) > 0 ){

    # compute vectors of step starting points with the refpoint
    b <- x[from[ft],-1,drop=FALSE]
    b <- t( apply( b, 1, function(x) x - p ) )

    # Compute length of this vector (= distance)
    blen <- sqrt(.rowSums(b^2, nrow(b),ncol(b)))
    r[ft] <- blen
  }
  r
}



#' Distance between pairs of tracks at every timepoint
#'
#' For every timepoint in the dataset, compute pairwise distances between coordinates.
#'
#' @param X a tracks object
#' @param searchRadius if specified, return only pairs that are within this distance of each other. 
#' Defaults to \code{Inf}, so if left unspecified, all pairs are returned.
#' @param times (optional) a vector of timePoints to check pairs at; by default this is just everything.
#' @param quietly (default FALSE) if TRUE, suppress warnings when there are no tracks with
#' overlapping timepoints and an empty dataframe is returned. 
#'
#' @return a dataframe with the following columns:
#' \describe{
#'   \item{cell1}{the id of the track to which the first coordinate belongs}
#'   \item{cell2}{the id of the track to which the second coordinate belongs}
#' 	 \item{t}{the time point at which their distance is assessed}
#' 	 \item{dist}{the distance between the coordinates at this time}
#' }
#'
#' @examples
#' ## compute find timepoints where two t cells are within 1 micron of each other.
#' pairsByTime( TCells, searchRadius = 1 )
#' 
#' ## indeed, the following two cells nearly touch:
#' plot( TCells[ c("24","9258") ] )
#' @export
pairsByTime <- function( X, searchRadius = Inf, times = timePoints(X), quietly = FALSE )
{
  # check X tracks object
  if( !is.tracks( X ) ) stop( "X must be a tracks object!" )
  if( length(X) == 0 ){
  	if( !quietly ) warning( "pairsByTime: tracks object X is empty. Returning an empty dataframe." )
  	return( data.frame() )
  }
  
  # check times; must be a subset of timePoints(X).
  if( !all( times %in% timePoints(X) ) ) stop( "times must all occur in X!" )

  # tracks to dataframe, split into a list with one dataframe per 
  # timepoint. In that dataframe, set the cellid as rownames and
  # keep only the coordinates.
  df <- as.data.frame( X )
  dataByTime <- split( df, df$t )
  
  # filter data for the timepoints in 'times' (by default, that's everything)
  dataByTime <- dataByTime[ as.character( times ) ]
  
  coordsByTime <- lapply( dataByTime, function(x){
    rownames(x) <- x$id
    keep.cols <- !is.element( colnames(x), c("id"))
    x <- x[,keep.cols]
    return(x)
  })
  # Per timepoint in the list, compute distance matrix of all the cells present
  # at that time, and then return pairwise distances in a dataframe.
  pbt <- lapply( coordsByTime, function(x){
    
    # distance from distance matrix
    # (get coordinates of all points at this time; compute distance matrix)
    coords <- x[,colnames(x)!="t"]# remove first column containing time.
    dm <- as.matrix( stats::dist( coords ))
    
    # rownames contain cell ids. Get all pairs with combn. This works only if 
    # there are at least two ids; return an empty dataframe otherwise.
    
    if( length( rownames(dm) ) > 1 ){
    	pairs <- as.data.frame( t(utils::combn( rownames(dm), 2 )) )
    	colnames(pairs) <- c( "cell1","cell2" )
    	
    	# add time info
		pairs$t <- x[1,"t"]
		pairs$dist <- apply( pairs, 1, function(y) dm[ y[1], y[2] ] )
	
		# filter pairs within distance searchRadius of each other
		pairs <- pairs[ pairs$dist <= searchRadius, ]
    } else {
    	pairs <- data.frame( cell1 = character(0), cell2 = character(0), t = numeric(0), dist = numeric(0) )
    }    
    return(pairs)
  })
  out <- do.call( rbind, pbt )
  if( nrow(out) > 0 ){
  	rownames(out) <- paste0( out[,1], "-", out[,2], "-", out[,3] )
  } else{
  	if( !quietly ) warning( "pairsByTime: no tracks share time points; returning an empty dataframe.")
  }
  return( out )
}

#' Angle between Two Steps
#'
#' Compute the angle between two steps in the dataset that occur at the same timepoint.
#'
#' @param X a tracks object
#' @param trackids a vector of two indices specifying the tracks to get steps from, or
#' a dataframe/matrix of two columns (where every row contains a pair of trackids to compute 
#' a step angle for)
#' @param t the timepoint at which the steps should start, or a vector of timepoints if trackids
#' is a matrix with multiple step pairs to compute angles for.
#' @param degrees logical; should angle be returned in degrees instead of radians? (defaults to \code{TRUE})
#' @param quietly logical; should a warning be returned if one or both of the steps are missing
#' in the data and the function returns NA?
#'
#' @return A single angle, or NA if the desired step is missing for one or both
#' of the tracks. If trackids is a matrix with multiple step pairs to compute angles for,
#' the output is a numeric vector of angles (or NA values).
#'
#' @seealso \code{\link{distanceSteps}} to compute the distance between the step starting
#' points, \code{\link{timePoints}} to list all timepoints in a dataset,
#' and \code{\link{AngleAnalysis}} for other methods to compute angles and distances.
#'
#' @examples
#' ## Find the angle between the steps of the tracks with ids 1 and 2, at the 3rd
#' ## timepoint in the dataset.
#' t <- timePoints( TCells )[3]
#' angleSteps( TCells, c("1","3"), t )
#' 
#' ## Do this for multiple pairs and times at once: between cells 1 and 3 at the
#' ## 3rd timepoint, and between 1 and 4 at the fourth timepoint.
#' pairs <- data.frame( cell1 = c("1","1"), cell2 = c("3","4"))
#' times <- timePoints(TCells)[3:4]
#' angleSteps( TCells, pairs, times )
#' @export
angleSteps <- function( X, trackids, t, degrees = TRUE, quietly = FALSE )
{
  # ---- Checking/handling input
  if( !is.tracks(X) ){
    stop( "angleSteps: X must be a tracks object." )
  }
  if( any( !is.element( t, timePoints(X) ) ) ){
    stop( "angleSteps: the supplied timepoints must correspond to those in the data!")
  }
  if( is.data.frame(trackids)) trackids <- as.matrix( trackids )
  if( is.matrix( trackids ) ){
    if( ncol(trackids) != 2 ){
      stop( "angleSteps: an angle is only defined for exactly 2 steps. Please provide exactly 2 trackids, or a matrix with 2 trackid columns per row.")
    }
    if( nrow( trackids ) != length(t) ){
      stop( "angleSteps: if angles are to be computed for a matrix of trackid pairs, then t must be a vector with the same length as the number of pairs.")
    }
  } else {
    if( !length(trackids) == 2 ){
      stop( "angleSteps: an angle is only defined for exactly 2 steps. Please provide exactly 2 trackids.")
    }
  }
  if( any( !is.element( trackids, names(X) ) ) ){
    stop( "angleSteps: cannot find all supplied trackids in the data.")
  }
  
  # if trackids is not a matrix (so we're checking just one angle), 
  # then just make it a matrix of 1 row.
  if( !is.matrix(trackids ) )  trackids <- matrix( trackids, ncol = length( trackids ) )
  
  # ---- Compute angle(s)
  # To compute angle(s), first restructure tracks to a dataframe; this will be faster
  # because it will allow to compute multiple angles at once using matrix-based computations.
  df <- as.data.frame(X)
  
  # Add a column with the index of the timepoint (rather than absolute time);
  # used below to find 'steps' in the data from coordinates at subsequent timepoints.
  tsort <- sort( timePoints(X) )
  tIndices <- stats::setNames( seq_along( tsort ), tsort )
  df$timeIndex <- tIndices[ as.character( df$t  ) ]
  
  # set rownames of the data to a tag containing the cellid and the timepoint; this will
  # allow us to find coordinates rapidly below.
  rownames(df) <- paste0( df$id, "-", df$timeIndex )
  
  # convert the timepoints of interest to time indices.
  tindexVec <- tIndices[ as.character(t) ]
  
  # Now find all relevant coordinates of the step pairs at the indicated time point.
  # In trackids, each row is a pair with the two trackids in the columns. For both cells
  # (columns), first find its coordinates at time t (="start") by using the rownames in df; 
  # discard t/id columns to keep coordinates only. 
  coordCols <- !is.element( colnames(df), c("t","id","timeIndex"))
  start1 <- df[ paste0( trackids[,1], "-", unname( tindexVec ) ), coordCols ]
  start2 <- df[ paste0( trackids[,2], "-", unname( tindexVec ) ), coordCols ]
  
  # also find the coordinates at the next timepoint ( = "end" of the step). Note that this
  # is simply the point with a time-index incremented by 1:
  end1 <- df[ paste0( trackids[,1], "-", unname( tindexVec+1 ) ), coordCols ]
  end2 <- df[ paste0( trackids[,2], "-", unname( tindexVec+1 ) ), coordCols ]
  
  # check for which ids all relevant coordinates are actually found
  keep <- ( rowSums( is.na( start1 )  ) == 0 ) & ( rowSums( is.na( start2 )  ) == 0 ) & ( rowSums( is.na( end1 )  ) == 0 ) & ( rowSums( is.na( end2 )  )== 0 )
  
  # subtract start from end to get displacement vector
  d1 <- as.matrix( end1[keep,] - start1[keep,] )
  d2 <- as.matrix( end2[keep,] - start2[keep,] )
  
  # init angles with NA (the default for if steps not in data), then add angles where possible
  ang <- rep( NA, nrow( trackids ) ) 
  if( sum( keep ) > 0 ) ang[keep] <- vecAngle( d1, d2, degrees = degrees )
  
  # warn if NAs left
  if( any( is.na( ang ) ) && !quietly ){
    warning( "Warning: for some pairs I cannot find data for both steps at the indicated time. Returning NA.")
  }
  
  ang
}

#' Distance between Two Steps
#'
#' Compute the distance between two steps in the dataset that occur at the same timepoint.
#' The distance is the distance between the step starting points.
#'
#' @param X a tracks object
#' @param trackids a vector of two indices specifying the tracks to get steps from, or
#' a dataframe/matrix of two columns (where every row contains a pair of trackids to compute 
#' a step angle for)
#' @param t the timepoint at which the steps should start, or a vector of such timepoints if
#' multiple step pairs are supplied in \code{trackids}.
#' @param quietly logical; should a warning be returned if one or both of the steps are missing
#' in the data and the function returns NA?
#'
#' @return A single distance (NA if the desired timepoint is missing for one or both
#' of the tracks), or a vector of such distances if multiple step pairs are supplied in \code{trackids}.
#'
#' @seealso \code{\link{angleSteps}} to compute the angle between the steps,
#' \code{\link{timePoints}} to list all timepoints in a dataset,
#' and \code{\link{AngleAnalysis}} for other methods to compute angles and distances.
#'
#' @examples
#' ## Find the distance between the steps of the tracks with ids 1 and 3, at the 3rd
#' ## timepoint in the dataset.
#' t <- timePoints( TCells )[3]
#' distanceSteps( TCells, c("1","3"), t )
#' 
#' ## Do this for multiple pairs and times at once: between cells 1 and 3 at the
#' ## 3rd timepoint, and between 1 and 4 at the fourth timepoint.
#' pairs <- data.frame( cell1 = c("1","1"), cell2 = c("3","4"))
#' times <- timePoints(TCells)[3:4]
#' distanceSteps( TCells, pairs, times )
#' @export
distanceSteps <- function( X, trackids, t, quietly = FALSE )
{
  
  # ---- Checking/handling input
  if( !is.tracks(X) ){
    stop( "distanceSteps: X must be a tracks object." )
  }
  if( any( !is.element( t, timePoints(X) ) ) ){
    stop( "distanceSteps: the supplied timepoints must correspond to those in the data!")
  }
  if( is.data.frame(trackids)) trackids <- as.matrix( trackids )
  if( is.matrix( trackids ) ){
    if( ncol(trackids) != 2 ){
      stop( "distanceSteps: a distance is only defined for exactly 2 steps. Please provide exactly 2 trackids, or a matrix with 2 trackid columns per row.")
    }
    if( nrow( trackids ) != length(t) ){
      stop( "distanceSteps: if distances are to be computed for a matrix of trackid pairs, then t must be a vector with the same length as the number of pairs.")
    }
  } else {
    if( !length(trackids) == 2 ){
      stop( "distanceSteps: only defined for exactly 2 steps. Please provide exactly 2 trackids.")
    }
  }
  if( any( !is.element( trackids, names(X) ) ) ){
    stop( "distanceSteps: cannot find all supplied trackids in the data.")
  }
  
  # if trackids is not a matrix (so we're checking just one angle), 
  # then just make it a matrix of 1 row.
  if( !is.matrix(trackids ) )  trackids <- matrix( trackids, ncol = length( trackids ) )
  
  # ---- Compute distance(s)
  # Compute end times of the steps
  tsort <- sort( timePoints(X) )
  tIndices <- stats::setNames( seq_along( tsort ), tsort )
  tIndexVec <- tIndices[as.character(t)]
  tNext <- tsort[ tIndexVec + 1 ]
  
  # check for which pair/time combinations step data are actually present
  df <- as.data.frame( X )
  rownames(df) <- paste0( df$id, "-", df$t )
  coordCols <- !is.element( colnames(df), c("t","id") )
  start1 <- df[ paste0( trackids[,1], "-", t ), coordCols ]
  start2 <- df[ paste0( trackids[,2], "-", t ), coordCols ]
  end1 <- df[ paste0( trackids[,1], "-", tNext ), coordCols ]
  end2 <- df[ paste0( trackids[,2], "-", tNext ), coordCols ]
  keep <- ( rowSums( is.na( start1 )  ) == 0 ) & ( rowSums( is.na( start2 )  ) == 0 ) & ( rowSums( is.na( end1 )  ) == 0 ) & ( rowSums( is.na( end2 )  )== 0 )
  
  # for the pairs where this is the case, we can get the distance using pairsByTime, using
  # the starting time of the step. For other pairs, we'll return NA, so init this first.
  distances <- rep( NA, nrow( trackids ) )
  if( sum( keep ) > 0 ){
    tvec <- intersect( unique( c( t, tNext ) ), timePoints( X ) )
    pbt <- pairsByTime( X, times = tvec )
    #trackids must be sorted with lowest first
    trackids <- t( apply( trackids, 1, sort ) ) 
    distances[keep] <- pbt[ paste0( trackids[keep,1], "-", trackids[keep,2], "-", t[keep] ),"dist" ]
  }
  
  # if trackids are equal their distance is zero.
  same <- ( trackids[,1] == trackids[,2] )
  distances[same] <- 0
  
  
  if( any(is.na(distances)) ){
    if(!quietly){warning( "Warning: cannot find data for both steps. Returning NA.")}
  }
  return(distances)
  
  
}


#' Find Pairs of Steps Occurring at the Same Time
#'
#' Find cell indices and timepoints where these cells both have a step.
#'
#' @param X a tracks object
#' @param filter.steps optional: a function used to filter steps on. See examples.
#'
#' @return A dataframe with three columns: two for the indices of cellpairs that
#' share a step, and one for the timepoint at which they do so.
#'
#' @examples
#' ## Find all pairs of steps in the T cell data that displace at least 2 microns.
#' pairs <- stepPairs( TCells, filter.steps = function(t) displacement(t) > 2 )
#'
#' @export
stepPairs <- function( X, filter.steps=NULL )
{

  if( !is.tracks(X) ){
    stop( "stepPairs: X must be a tracks object!" )
  }

  dout <- data.frame()
  for( t in timePoints(X) ){
    # All steps starting at that timepoint
    xt <- subtracksByTime( X, t, 1 )

    # Filter out steps that do not fulfill the filter criterion
    if( length(xt) > 0 && is.function( filter.steps ) ){
      xt <- xt[ sapply(xt,filter.steps ) ]
    }

    # Extract their ids
    ids <- names(xt)

    # Make all possible pairs
    if( length(ids) >= 2 ){
      pairs <- as.data.frame( t( utils::combn( ids, 2 ) ),
                              stringsAsFactors = FALSE )
      colnames(pairs) <- c( "p1","p2" )

      # To a dataframe
      pairs$t <- t
      dout <- rbind( dout, pairs )
    }
  }
  dout
}



#' Find Distances and Angles for all Pairs of Steps
#'
#' Find cell indices and timepoints where these cells both have a step, then return
#' angles and distances for each pair of steps.
#'
#' @param X a tracks object
#' @param searchRadius if specified, only return analysis for pairs of steps that
#' start within distance searchRadius from each other
#' @param filter.steps optional: a function used to filter steps on. See examples.
#' @param quietly (default FALSE) if TRUE, suppress warnings
#' @param ... further arguments passed on to \code{angleSteps}
#'
#' @return A dataframe with five columns: two for the indices of cellpairs that
#' share a step, one for the timepoint at which they do so, one for the distance
#' between them, and one for their angle.
#'
#' @details Analyzing step angles at different distances can be useful to detect
#' directional bias or local crowding effects; see (Beltman et al, 2009).
#'
#' Internally, the function uses \code{\link{stepPairs}}, \code{\link{angleSteps}},
#' and \code{\link{distanceSteps}}.
#'
#' @seealso \code{\link{analyzeCellPairs}} to do something similar for entire tracks
#' rather than single steps.
#'
#'
#' @references
#' Joost B. Beltman, Athanasius F.M. Maree and Rob. J. de Boer (2009),
#' Analysing immune cell migration. \emph{Nature Reviews Immunology} \bold{9},
#' 789--798. doi:10.1038/nri2638
#'
#' @examples
#' ## Plot distance versus angle for all step pairs, filtering for those that
#' ## displace at least 2 microns. Sample dataset in this example for speed.
#' pairs <- analyzeStepPairs( sample( TCells, 100), filter.steps = function(t) displacement(t) > 2 )
#' scatter.smooth( pairs$dist, pairs$angle )
#' @export
analyzeStepPairs <- function( X, filter.steps = NULL, searchRadius = Inf, quietly = FALSE, ... )
{
  
  if( !is.tracks(X) ){
    stop( "analyzeStepPairs: X must be a tracks object!" )
  }
  
  # Obtain cell pairs for each timepoint; using filter.steps as filter criterion.
  pairs <- stepPairs( X, filter.steps = filter.steps )
  if( nrow(pairs) == 0 ) return(pairs)
  
  # initialise dist/angle with NA. Set rownames p1-p2-t for lookup later.
  pairs$dist <- NA
  pairs$angle <- NA
  rownames(pairs) <- paste0( pairs$p1, "-", pairs$p2, "-", pairs$t )
  
  # pairs of cells by time coordinate; this also computes a distance. 
  # these data have rownames id1-id2-t for matching with the pairs dataframe later.
  pByTime <- pairsByTime( X, searchRadius = searchRadius, quietly = quietly )
  
  # add angles for these pairs (which are within the search radius).
  pByTime$angle <- angleSteps( X, as.matrix(pByTime[,1:2]), pByTime[,3], quietly = TRUE, ... )
  
  # Note that pByTime contains all pairs that co-occur at some timepoint, but that
  # for analyzing steps they need to co-occur at least at two subsequent timepoints.
  # Thus, pByTime may contain pairs not included in the 'pairs' dataframe generated 
  # by stepPairs() above. The same holds true if some filter.steps criterion is specified.
  # But it is also possible that there are pairs in 'pairs' which are not in pByTime,
  # since pByTime is already focused on pairs within distance searchRadius from each other.
  # therefore intersect to keep only the pairs that fulfill all these criteria. 
  pairs <- pairs[ intersect( rownames(pairs), rownames( pByTime ) ), ]
  
  # Now that we are sure all pairs in 'pairs' are also in 'pByTime', lookup the computed
  # angles and distances there and return. 
  pairs$dist <- pByTime[ rownames(pairs), "dist" ]
  pairs$angle <- pByTime[ rownames(pairs), "angle" ]
  pairs
}



#' Find Distances and Angles for all Pairs of Tracks
#'
#' Find all pairs of cells and return the shortest distance between them at any
#' point in time (if they share any time points), as well as the angle between 
#' their overall displacement vectors.
#'
#' @param X a tracks object
#' @param searchRadius if specified, only return analysis for pairs of cells that
#' are within distance searchRadius from each other at least at one point in time.
#' @param quietly (default FALSE) if TRUE, suppress warnings
#' @param ... further arguments passed on to \code{angleCells}
#'
#' @return A dataframe with four columns: two for the indices of cellpairs,
#' one for the distance between them, and one for their angle. Note that the 
#' distance will be NA for pairs of tracks that do not share time points, but
#' their angle will still be computed.
#'
#' @details Analyzing track angles at different distances can be useful to detect
#' directional bias or local crowding effects; see (Beltman et al, 2009).
#'
#' Internally, the function uses \code{\link{cellPairs}}, \code{\link{angleCells}},
#' and \code{\link{distanceCells}}.
#'
#' @seealso \code{\link{analyzeStepPairs}} to do something similar for single steps
#' rather than entire tracks.
#'
#'
#' @references
#' Joost B. Beltman, Athanasius F.M. Maree and Rob. J. de Boer (2009),
#' Analysing immune cell migration. \emph{Nature Reviews Immunology} \bold{9},
#' 789--798. doi:10.1038/nri2638
#'
#' @examples
#' ## Plot distance versus angle for all cell pairs. Sample T-cell data here for speed.
#' pairs <- analyzeCellPairs( sample( TCells, 100 ) )
#' scatter.smooth( pairs$dist, pairs$angle )
#' @export
analyzeCellPairs <- function( X, searchRadius = Inf, quietly = FALSE, ... )
{
  
  if( !is.tracks(X) ){
    stop( "analyzeCellPairs: X must be a tracks object!" )
  }
  
  # Make all possible pairs of cellids
  pairs <- cellPairs( X )
  if( nrow( pairs ) == 0 ) return(pairs)
  
  # init distances and angles with NA.
  pairs$dist <- NA
  pairs$angle <- NA
  
  # per timepoint, compute distances for all pairs of cells that
  # are both present at that time point. Output is a df with cell ids, time, and distance.
  # Note that multiple distances per pair are present if the cells co-occur at multiple
  # timepoints.
  allPairs <- pairsByTime( X, searchRadius = searchRadius, quietly = quietly )
  if( nrow( allPairs ) == 0 ){
  	return( data.frame( cell1 = character(0), cell2 = character(0), dist = numeric(0), angle = numeric(0) ))
  }
  
  # Split by pair to find the min distance between the cells over time; then merge again
  # into a single dataframe.
  distByPair <- split( allPairs, paste0( allPairs$cell1, "-", allPairs$cell2 ) )
  minDistByPair <- lapply( distByPair, function(x) {
    minD <- min( x$dist )
    x <- x[1,]
    x$dist <- minD
    return(x)
  })
  pairMinDist <- do.call( rbind, minDistByPair )
  
  # the pairMinDist already has rownames in the form of 'id1-id2' so data can be accessed easily.
  # do the same with the 'pairs' dataframe of query points generated above:
  rownames( pairs ) <- paste0( pairs$cell1, "-", pairs$cell2 )
  
  # For all the query points, get the distance and angle from the pairMinDist dataframe.
  # Some query points will not be in this dataframe (if cells do not co-occur at any time);
  # in that case, dist will remain NA. 
  pairs[ rownames( pairMinDist ), "dist" ] <- pairMinDist$dist 
  
  # add angles and return
  pairs$angle <- angleCells( X, pairs[,1:2], ...)
  return(pairs)
  
}


#' Find Pairs of Tracks
#'
#' Get all unique combinations of two track ids.
#'
#' @param X a tracks object
#'
#' @return A dataframe with two columns: one for each of the track ids in the pair.
#' Each row represents a pair.
#'
#' @examples
#' ## Find all pairs of cells in the T cell data
#' pairs <- cellPairs( TCells )
#' @export
cellPairs <- function( X )
{
  if( !is.tracks(X) ){
    stop( "cellPairs: X must be a tracks object!" )
  }

  cellids <- names( X )
  pairs <- data.frame()

  if( length(cellids) >= 2 ){
    # Make all possible pairs of cellids
    pairs <- as.data.frame( t( utils::combn( cellids, 2 ) ),
                            stringsAsFactors = FALSE )
    colnames(pairs) <- c( "cell1","cell2" )
  }
  return(pairs)
}

#' Angle between Two Tracks
#'
#' Compute the angle between the displacement vectors of two tracks in the dataset,
#' or of several such pairs at once.
#' Note that in contrast to \code{\link{distanceCells}}, this angle is computed even
#' when the two tracks do not share any time points.
#'
#' @param X a tracks object
#' @param cellids a vector of two indices specifying the tracks to get steps from, or
#' a dataframe/matrix of two columns (where every row contains a pair of cellids to compute 
#' an angle for)
#' @param degrees logical; should angle be returned in degrees instead of radians? (defaults to \code{TRUE})
#'
#' @return A single angle (if two cellids given), or a vector of angles (if multiple pairs of cellids are supplied).
#'
#' @seealso \code{\link{distanceCells}} to compute the minimum distance between the tracks,
#' and \code{\link{AngleAnalysis}} for other methods to compute angles and distances.
#'
#' @examples
#' ## Find the angle between the tracks with ids 1 and 3
#' angleCells( TCells, c("1","3") )
#' 
#' ## Find the angles of several cell pairs at once
#' pairs <- data.frame( cell1 = c("1","1"), cell2 = c( "3","4" ) )
#' angleCells( TCells, pairs )
#' @export
angleCells <- function (X, cellids, degrees = TRUE ) 
{
  #--- check/handle inputs
  if (!is.tracks(X)) {
    stop("angleCells: X must be a tracks object!")
  }
  if( is.data.frame(cellids)) cellids <- as.matrix( cellids )
  unique.ids <- unique( as.vector( cellids ))
  if (any(!is.element(unique.ids, names(X)))) {
    stop("angleCells: cannot find all cellids in data.")
  }
  if( is.matrix( cellids ) ){
    if( ncol(cellids) != 2 ){
      stop( "angleCells: an angle is only defined for exactly 2 cells. Please provide exactly 2 cellids, or a matrix with 2 trackid columns per row.")
    }
  } else {
    if( !length(cellids) == 2 ){
      stop( "angleCells: an angle is only defined for exactly 2 cells. Please provide exactly 2 cellids.")
    }
  }
  
  
  # if cellids is not a matrix (so we're checking just one angle), 
  # then just make it a matrix of 1 row.
  if( !is.matrix(cellids ) )  cellids <- matrix( cellids, ncol = length(cellids ) )
  
  #--- compute
  # Get displacement vectors for all the ids in trackids
  disps <- lapply( unique( as.vector( cellids ) ), function(x) displacementVector( X[[x]] )  )
  names(disps) <- unique( as.vector( cellids ) )
  
  # Get a matrix of displacement vectors for the ids in the first column of trackids, 
  # and likewise for the second id of each pair
  disps1 <- do.call( rbind, disps[ cellids[,1] ] )
  disps2 <- do.call( rbind, disps[ cellids[,2] ] )
  
  # compute all angles at once using these matrices.
  return(vecAngle(disps1, disps2, degrees = degrees))
}


#' Minimum Distance between Two Cells
#'
#' Compute the minimum distance between two cells in the dataset (minimum over all)
#' the timepoints where they were both measured.
#'
#' @param X a tracks object
#' @param cellids a vector of two indices specifying the tracks to compute distance between, or
#' a dataframe/matrix of two columns (where every row contains a pair of cellids to compute 
#' a distance for)
#' @param quietly if TRUE, suppress warnings about returning NA distances.
#' 
#' @return A single distance (NA if the the tracks do not have overlapping timepoints), or 
#' a vector of such distances if multiple pairs are supplied in \code{cellids}.
#'
#' @seealso \code{\link{angleCells}} to compute the angle between the track displacement vectors,
#' and \code{\link{AngleAnalysis}} for other methods to compute angles and distances.
#'
#' @examples
#' ## Find the minimum distance between the tracks with ids 1 and 3
#' distanceCells( TCells, c("1","3") )
#' @export
distanceCells <- function( X, cellids, quietly = FALSE )
{
  if( !is.tracks(X) ){
    stop( "distanceCells: X must be a tracks object!" )
  }
  if( is.data.frame(cellids)) cellids <- as.matrix( cellids )
  if( any( !is.element( unique( as.vector( cellids ) ) , names(X) ) ) ){
    stop( "distanceCells: cannot find both cellids in data." )
  }
  if( is.matrix( cellids ) ){
    if( ncol(cellids) != 2 ){
      stop( "distanceCells: a distance is only defined for exactly 2 cells. Please provide exactly 2 cellids, or a matrix with 2 trackid columns per row.")
    }
  } else {
    if( !length(cellids) == 2 ){
      stop( "distanceCells: a distance is only defined for exactly 2 cells. Please provide exactly 2 cellids.")
    }
  }
  
  # if cellids is not a matrix (so we're checking just one angle), 
  # then just make it a matrix of 1 row.
  if( !is.matrix(cellids ) )  cellids <- matrix( cellids, ncol = length(cellids ) )
  
  # Sort cellids per row
  cellids <- t( apply( cellids, 1, sort ) )


  # initalize distance as NA. 
  pairs <- paste0( cellids[,1],"-", cellids[,2] )
  distances <- stats::setNames( rep( NA, length(pairs)), pairs )

  # for any pairs that share timepoints, add the distance. 
  # pairsByTime returns distances for all timepoints where both tracks are present;
  # use this to extract the min distance for each pair. 
  pbt <- pairsByTime( X[ unique( as.vector( cellids ) ) ], quietly = TRUE )
  if( nrow(pbt) > 0 ){
  	pbt[,1:2] <- t( apply( pbt[,1:2], 1, sort ) ) 
  	pbtMin <- sapply( split( pbt, paste0( pbt$cell1, "-", pbt$cell2 ) ), function(x) min( x$dist ) )
  	add.pairs <- pairs[ is.element( pairs, names( pbtMin ))]
  	distances[add.pairs] <- unname( pbtMin[add.pairs] )
  }
  
  # if cellids are equal, distance is zero
  distances[ cellids[,1]==cellids[,2]] <- 0

  # Warning if there are NAs
  if( any(is.na(distances)) ){
    if(!quietly){warning( "Warning: distance undefined for cells that don't share timepoints; returning NA.")}
  }
  return(unname( distances) )

}

