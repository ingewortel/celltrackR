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
#'  reference point. The angle is between the direction of the first step in the track and the
#'  line between its first coordinate and the reference point. Angles are by default returned in
#'  degrees, use \code{degrees=FALSE} to obtain radians. These functions are useful to detect
#'  directional bias towards a point of interest, which would result in an average angle of less
#'  than 90 degrees with the reference point (especially for tracks at a small distance to the
#'  reference point).
#'
#' \code{\link{angleToPlane}} and \code{\link{distanceToPlane}} return the angle/distance of the track to a
#'  plane of interest. This plane must be specified by three points lying on it.
#'  The distance returned is between the first coordinate in the track and the
#'  reference point. The angle is between the direction of the first step in the track and the
#'  plane of interest. These functions are useful to detect tracking artefacts near the borders
#'  of the imaging volume. Use \code{\link{boundingBox}} to guess where those borders are.
#'  Angles are by default returned in
#'  degrees, use \code{degrees=FALSE} to obtain radians.
#'
#' \code{\link{angleToDir}} returns the angle of the first step in a track to a direction of interest.
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
#' steps <- subtracks( Neutrophils, 1 )
#' bb <- boundingBox( Neutrophils )
#' angles <- sapply( steps, angleToPoint, p = bb["max",-1] )
#' distances <- sapply( steps, distanceToPoint, p = bb["max",-1] )
#' scatter.smooth( distances, angles )
#' abline( h = 90, col = "red" )
#'
#' ## Get a distribution of Neutrophil step angles with the reference direction
#' ## in positive y direction. The histogram is enriched for low angles, suggesting
#' ## directed movement:
#' hist( sapply( steps, angleToDir, dvec=c(0,1,0) ) )
#'
#' ## Plotting the angle versus the distance to a reference plane can be informative to
#' ## detect tracking artefacts near the border of the imaging volume.
#' ## We should be suspicious especially when small angles are more frequent at low distances
#' ## to the border planes.
#' steps <- subtracks( TCells, 1 )
#' minz <- boundingBox( TCells )["min","z"]
#' ## Compute angles and distances to the lower plane in z-dimension
#' angles <- sapply( steps, angleToPlane, p1 = c(0,0,minz), p2 = c(1,0,minz), p3 = c(0,1,minz) )
#' distances <- sapply( steps, distanceToPlane, p1 = c(0,0,minz), p2 = c(1,0,minz), p3 = c(0,1,minz) )
#' scatter.smooth( distances, angles )
#' abline( h = 32.7, col = "red" )
#'
#' ## Plot distance versus angle for all cell pairs
#' pairs <- analyzeCellPairs( TCells )
#' scatter.smooth( pairs$dist, pairs$angle )
#' abline( h = 90, col = "red" )
#'
#' ## Plot distance versus angle for all step pairs, filtering for those that
#' ## displace at least 2 microns
#' pairs <- analyzeStepPairs( TCells, filter.steps = function(t) displacement(t) > 2 )
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
#' Compute the angle between the first step of a track and a reference point. Useful to
#' detect directed movement towards a point (see examples).
#'
#' @param x a single input track; a matrix whose first column is time and whose
#'  remaining columns are a spatial coordinate.
#' @param p numeric vector of coordinates of the reference point p to compute angles/distances to.
#' @param from index, or vector of indices, of the first row of the track. If
#' \code{from} is a vector, angles are returned for all steps starting at
#' the indices in \code{from}.
#' @param xdiff row differences of x.
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
#' ## Use bb to get the corner with highest x,y,and z value
#' ## The histogram is enriched for low angles, suggesting directed movement:
#' steps <- subtracks( Neutrophils, 1 )
#' bb <- boundingBox( Neutrophils )
#' hist( sapply( steps, angleToPoint, p = bb["max",-1] ) )
#'
#' ## The same does not hold for movement of T cells towards the point (0,0,0)
#' steps <- subtracks( TCells, 1 )
#' hist( sapply( steps, angleToPoint, p = c(0,0,0) ) )
#'
#' ## Plotting the angle versus the distance to the reference point can also be informative,
#' ## especially when small angles are more frequent at lower distances.
#' angles <- sapply( steps, angleToPoint, p = bb["max",-1] )
#' distances <- sapply( steps, distanceToPoint, p = bb["max",-1] )
#' scatter.smooth( distances, angles )
#' abline( h = 90, col = "red" )
#' @export
angleToPoint <- function (x, p = c(1,1,1), from = 1, xdiff = diff(x), degrees = TRUE )
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
    a <- xdiff[from[ft], -1, drop = FALSE]

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
#' Compute the angle between the first step of a track and a reference direction.
#' Useful to detect biased movement when the directional bias is known (see examples).
#'
#' @param x a single input track; a matrix whose first column is time and whose
#'  remaining columns are a spatial coordinate.
#' @param dvec numeric vector specifying a reference direction to compute angles to.
#' @param from index, or vector of indices, of the first row of the track. If
#' \code{from} is a vector, angles are returned for all steps starting at
#' the indices in \code{from}.
#' @param xdiff row differences of x.
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
#' hist( sapply( steps, angleToDir, dvec=c(0,1,0) ) )
#' @export
angleToDir <- function (x, dvec = c(1,1,1), from = 1, xdiff = diff(x), degrees=TRUE )
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
    a <- xdiff[from[ft], -1, drop = FALSE]
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
#' Compute the angle between the first step of a track and a reference plane.
#' Useful to detect directed movement and/or tracking artefacts.
#'
#' @param x a single input track; a matrix whose first column is time and whose
#'  remaining columns are a spatial coordinate.
#' @param p1,p2,p3 numeric vectors of coordinates of three points specifying a reference plane to
#'  compute distances to.
#' @param from index, or vector of indices, of the first row of the track. If
#' \code{from} is a vector, angles are returned for all steps starting at
#' the indices in \code{from}.
#' @param xdiff row differences of x.
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
#' steps <- subtracks( TCells, 1 )
#' minz <- boundingBox( TCells )["min","z"]
#' ## Compute angles and distances to the lower plane in z-dimension
#' angles <- sapply( steps, angleToPlane, p1 = c(0,0,minz), p2 = c(1,0,minz), p3 = c(0,1,minz) )
#' distances <- sapply( steps, distanceToPlane, p1 = c(0,0,minz), p2 = c(1,0,minz), p3 = c(0,1,minz) )
#' scatter.smooth( distances, angles )
#' abline( h = 32.7, col = "red" )
#' @export
angleToPlane <- function (x, p1 = c(0,0,0), p2 = c(0,1,0), p3 = c(1,0,0),
                          from = 1, xdiff = diff(x), degrees =TRUE )
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
    a <- xdiff[from[ft], -1, drop = FALSE]

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
#' steps <- subtracks( TCells, 1 )
#' minz <- boundingBox( TCells )["min","z"]
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

#' Angle between Two Steps
#'
#' Compute the angle between two steps in the dataset that occur at the same timepoint.
#'
#' @param X a tracks object
#' @param trackids a vector of two indices specifying the tracks to get steps from.
#' @param t the timepoint at which the steps should start.
#' @param degrees logical; should angle be returned in degrees instead of radians? (defaults to \code{TRUE})
#' @param quietly logical; should a warning be returned if one or both of the steps are missing
#' in the data and the function returns NA?
#'
#' @return A single angle, or NA if the desired timepoint is missing for one or both
#' of the tracks.
#'
#' @seealso \code{\link{distanceSteps}} to compute the distance between the step starting
#' points, \code{\link{timePoints}} to list all timepoints in a dataset,
#' and \code{\link{AngleAnalysis}} for other methods to compute angles and distances.
#'
#' @examples
#' ## Find the angle between the steps of the tracks with ids 1 and 2, at the 3rd
#' ## timepoint in the dataset.
#' t <- timePoints( TCells )[3]
#' angleSteps( TCells, c("1","2"), t )
#' @export
angleSteps <- function( X, trackids, t, degrees = TRUE, quietly = FALSE )
{
  if( !is.tracks(X) ){
    stop( "angleSteps: X must be a tracks object." )
  }
  if( !length(trackids) == 2 ){
    stop( "angleSteps: an angle is only defined for exactly 2 steps. Please provide exactly 2 trackids.")
  }
  if( any( !is.element( trackids, names(X) ) ) ){
    stop( "angleSteps: cannot find all supplied trackids in the data.")
  }
  X2 <- selectSteps( X, trackids, t )
  if( any( sapply(X2, is.null ) ) ){
    if(!quietly){warning( "Warning: cannot find data for both steps. Returning NA.")}
    return(NA)
  }
  if( length(X2) != 2 ){
    if(!quietly){warning( "Warning: cannot find data for both steps. Returning NA.")}
    return(NA)
  }

  a <- diff( X2[[1]][,-1] )
  b <- diff( X2[[2]][,-1] )

  # # Normalize a and b by their length
  # a <- a/sqrt(sum(a^2))
  # b <- b/sqrt(sum(b^2))
  #
  # # Dot product of normalized a and b; this equals cos alpha
  # rs <- sum(a * b)
  #
  # # angle is acos.
  # ang <- acos(rs)
  ang <- vecAngle( a, b, degrees = degrees )
  ang
}

#' Distance between Two Steps
#'
#' Compute the distance between two steps in the dataset that occur at the same timepoint.
#' The distance is the distance between the step starting points.
#'
#' @param X a tracks object
#' @param trackids a vector of two indices specifying the tracks to get steps from.
#' @param t the timepoint at which the steps should start.
#' @param quietly logical; should a warning be returned if one or both of the steps are missing
#' in the data and the function returns NA?
#'
#' @return A single distance, or NA if the desired timepoint is missing for one or both
#' of the tracks.
#'
#' @seealso \code{\link{angleSteps}} to compute the angle between the steps,
#' \code{\link{timePoints}} to list all timepoints in a dataset,
#' and \code{\link{AngleAnalysis}} for other methods to compute angles and distances.
#'
#' @examples
#' ## Find the distance between the steps of the tracks with ids 1 and 2, at the 3rd
#' ## timepoint in the dataset.
#' t <- timePoints( TCells )[3]
#' distanceSteps( TCells, c("1","2"), t )
#' @export
distanceSteps <- function( X, trackids, t, quietly = FALSE )
{

  if( !is.tracks(X) ){
    stop( "distanceSteps: X must be a tracks object." )
  }
  if( !length(trackids) == 2 ){
    stop( "distanceSteps: only defined for exactly 2 steps. Please provide exactly 2 trackids.")
  }
  if( any( !is.element( trackids, names(X) ) ) ){
    stop( "distanceSteps: cannot find all supplied trackids in the data.")
  }
  X2 <- selectSteps( X, trackids, t )
  if( any( sapply(X2, is.null ) ) ){
    if(!quietly){warning( "Warning: cannot find data for both steps. Returning NA.")}
    return(NA)
  }
  if( length(X2) != 2 ){
    if(!quietly){warning( "Warning: cannot find data for both steps. Returning NA.")}
    return(NA)
  }

  # Select the relevant timepoint, and extract coordinates (remove id/time columns)
  coords <- as.data.frame.tracks( subtracksByTime( X2, t, 0 ) )[,-c(1,2)]

  # Compute the euclidian distance
  diffm <- diff( as.matrix(coords) )
  return( sqrt( sum( diffm^2 ) ) )

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
#' @param filter.steps optional: a function used to filter steps on. See examples.
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
#' ## displace at least 2 microns
#' pairs <- analyzeStepPairs( TCells, filter.steps = function(t) displacement(t) > 2 )
#' scatter.smooth( pairs$dist, pairs$angle )
#' @export
analyzeStepPairs <- function( X, filter.steps = NULL, ... )
{

  if( !is.tracks(X) ){
    stop( "analyzeStepPairs: X must be a tracks object!" )
  }

  # Obtain cell paris for each timepoint
  pairs <- stepPairs( X, filter.steps = filter.steps )

  if( nrow(pairs) > 0 ){
    # Find the distance between the step starting points
    distances <- unname( apply( pairs, 1, function(x) distanceSteps( X, x[1:2], as.numeric(x[3]), quietly = TRUE ) ) )

    # Find the angles between the steps
    angles <- unname( apply( pairs, 1, function(x) angleSteps( X, x[1:2], as.numeric(x[3]), quietly = TRUE, ... ) ) )

    # Add to dataframe and return
    pairs$dist <- distances
    pairs$angle <- angles
  }
  pairs
}



#' Find Distances and Angles for all Pairs of Tracks
#'
#' Find all pairs of cells and return the shortest distance between them at any
#' point, as well as the angle between their overall displacement vectors.
#'
#' @param X a tracks object
#' @param ... further arguments passed on to \code{angleCells}
#'
#' @return A dataframe with four columns: two for the indices of cellpairs,
#' one for the distance between them, and one for their angle.
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
#' ## Plot distance versus angle for all cell pairs
#' pairs <- analyzeCellPairs( TCells )
#' scatter.smooth( pairs$dist, pairs$angle )
#' @export
analyzeCellPairs <- function( X, ... )
{

  if( !is.tracks(X) ){
    stop( "analyzeCellPairs: X must be a tracks object!" )
  }

  # Make all possible pairs of cellids
  pairs <- cellPairs( X )

  if( nrow( pairs ) > 0 ){
    # Compute angles and distances for all cell pairs in the data
    cellangles <- apply( pairs, 1, function(x)
      angleCells( X, x, ... ) )

    celldistances <- apply( pairs, 1, function(x)
      distanceCells( X, x ) )

    # Make a dataframe
    pairs$dist <- celldistances
    pairs$angle <- cellangles
  }
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
#' Compute the angle between the displacement vectors of two tracks in the dataset.
#'
#' @param X a tracks object
#' @param cellids a vector of two indices specifying the tracks to get steps from.
#' @param degrees logical; should angle be returned in degrees instead of radians? (defaults to \code{TRUE})
#'
#' @return A single angle.
#'
#' @seealso \code{\link{distanceCells}} to compute the minimum distance between the tracks,
#' and \code{\link{AngleAnalysis}} for other methods to compute angles and distances.
#'
#' @examples
#' ## Find the angle between the tracks with ids 1 and 2
#' angleCells( TCells, c("1","2") )
#' @export
angleCells <- function( X, cellids, degrees = TRUE )
{
  if( !is.tracks(X) ){
    stop( "angleCells: X must be a tracks object!" )
  }
  if( any( !is.element( cellids, names(X) ) ) ){
    stop( "angleCells: cannot find both cellids in data." )
  }
  X <- X[cellids]

  a <- displacementVector( X[[1]] )
  b <- displacementVector( X[[2]] )
  return( vecAngle( a, b, degrees = degrees ) )

}

#' Minimum Distance between Two Cells
#'
#' Compute the minimum distance between two cells in the dataset (minimum over all)
#' the timepoints where they were both measured.
#'
#' @param X a tracks object
#' @param cellids a vector of two indices specifying the tracks to compute distance between.
#'
#' @return A single distance, or NA if the the tracks do not have overlapping timepoints.
#'
#' @seealso \code{\link{angleCells}} to compute the angle between the track displacement vectors,
#' and \code{\link{AngleAnalysis}} for other methods to compute angles and distances.
#'
#' @examples
#' ## Find the minimum distance between the tracks with ids 1 and 2
#' distanceCells( TCells, c("1","2") )
#' @export
distanceCells <- function( X, cellids )
{
  if( !is.tracks(X) ){
    stop( "distanceCells: X must be a tracks object!" )
  }
  if( any( !is.element( cellids, names(X) ) ) ){
    stop( "distanceCells: cannot find both cellids in data." )
  }

  # Find timepoints occurring in both tracks
  t1 <- timePoints( X[ cellids[1] ] )
  t2 <- timePoints( X[ cellids[2] ] )
  t <- c( t1, t2 )
  tboth <- t[ duplicated(t) ]

  # If cells have no overlapping timepts, return NA.
  if( length(tboth) == 0 ){
    return(NA)
  }

  # Get the coordinate matrix of both cells for the overlapping timeinterval
  m1 <- X[[ cellids[1] ]]
  rownames( m1 ) <- m1[,"t"]
  m1 <- m1[ as.character(tboth), , drop =FALSE ]
  m2 <- X[[ cellids[2] ]]
  rownames( m2 ) <- m2[,"t"]
  m2 <- m2[ as.character(tboth), , drop = FALSE ]

  # Remove the time column from each matrix
  m1 <- m1[,-1, drop = FALSE]
  m2 <- m2[,-1, drop = FALSE]

  # Subtract from each other to get dx,dy,dz at each t
  mdiff <- m1 - m2

  # compute distance at each timepoint
  mdiff <- mdiff^2
  distances <- sqrt( rowSums( mdiff ) )
  return( min(distances) )


}


