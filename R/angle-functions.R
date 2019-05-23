#' Angle analysis
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
#'
#'
#' @name AngleAnalysis
#'
#' @seealso \code{\link{TrackMeasures}} for other measures that can be used to quantify tracks.
#'
#' See the vignettes on Quality Control and Track Analysis for more detailed examples of
#' angle analyses.
#' \code{browseVignettes( package = "celltrackR" )}
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
#' ## Get a distribution of Neutrophil step angles with the reference direction in positive y direction.
#' ## The histogram is enriched for low angles, suggesting directed movement:
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
vecAngle <- function( a, b, degrees = TRUE )
{
  # Check if inputs have the right dimensions
  if( class(a) != class(b) ){
    stop( "vecAngle: a and b must both be numeric vectors (of the same length) or matrices (of the same size)")
  }
  if( any( pracma::size(a) != pracma::size(b) ) ){
    stop( "vecAngle: cannot compute angle between vectors of unequal dimensions.")
  }
  if( is.matrix(a) && nrow(a) != nrow(b) ){
    stop( "vecAngle: a and b must have an equal number of rows.")
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

#' Angle With A Reference Point
#'
#' Compute the angle between the first step of a track and a reference point. Useful to
#' detect directed movement towards a point (see examples).
#'
#' @param x a single input track; a matrix whose first column is time and whose
#'  remaining columns are a spatial coordinate.
#' @param from index, or vector of indices, of the first row of the track. If
#' \code{from} is a vector, angles are returned for all steps starting at
#' the indices in \code{from}.
#' @param xdiff row differences of x.
#' @param degrees logical; should angles be returned in degrees rather than radians? (default = TRUE).
#' @param p numeric vector of coordinates of the reference point p to compute angles/distances to.
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
angleToPoint <- function (x, from = 1, p = c(1,1,1), xdiff = diff(x), degrees = TRUE )
{

  # Check if the given point has the correct dimensions
  if( length(p) != ncol(x) - 1 ){
    stop("In angleToPoint: Reference point coordinates must have the same number of dimensions as
         coordinates in the tracking data.")
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


#' Angle With A Reference Direction
#'
#' Compute the angle between the first step of a track and a reference direction.
#' Useful to detect biased movement when the directional bias is known (see examples).
#'
#' @param x a single input track; a matrix whose first column is time and whose
#'  remaining columns are a spatial coordinate.
#' @param from index, or vector of indices, of the first row of the track. If
#' \code{from} is a vector, angles are returned for all steps starting at
#' the indices in \code{from}.
#' @param xdiff row differences of x.
#' @param degrees logical; should angles be returned in degrees rather than radians? (default = TRUE).
#' @param dvec numeric vector specifying a reference direction to compute angles to.
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
#' ## Get a distribution of Neutrophil step angles with the reference direction in positive y direction.
#' ## The histogram is enriched for low angles, suggesting directed movement:
#' steps <- subtracks( Neutrophils, 1 )
#' hist( sapply( steps, angleToDir, dvec=c(0,1,0) ) )
angleToDir <- function (x, from = 1, dvec = c(1,1,1), xdiff = diff(x), degrees=TRUE )
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


#' Angle With A Reference Plane
#'
#' Compute the angle between the first step of a track and a reference plane.
#' Useful to detect directed movement and/or tracking artefacts.
#'
#' @param x a single input track; a matrix whose first column is time and whose
#'  remaining columns are a spatial coordinate.
#' @param from index, or vector of indices, of the first row of the track. If
#' \code{from} is a vector, angles are returned for all steps starting at
#' the indices in \code{from}.
#' @param xdiff row differences of x.
#' @param degrees logical; should angles be returned in degrees rather than radians? (default = TRUE).
#' @param p1,p2,p3 numeric vectors of coordinates of three points specifying a reference plane to
#'  compute distances to.
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
angleToPlane <- function (x, from = 1, p1 = c(0,0,0), p2 = c(0,1,0), p3 = c(1,0,0),
                          xdiff = diff(x), degrees =TRUE )
{

  # Check if the given points have the correct dimensions
  plengths <- c( length(p1), length(p2), length(p3) )
  if( length(unique(plengths)) > 1 ){
    stop("In angleToPlane: Points p1,p2,p3 specifying the plane must have the
         same number of coordinates.")
  }
  if( length(p1) != ( ncol(x) - 1 ) ){
    stop("In angleToPlane: Plane points must have the same number of coordinates as
         coordinates in the tracking data.")
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
    r[ft] <- vecAngle( a, pnorm, degrees = TRUE )

    # angle to the plane itself is (pi/2)-r
    deg90 <- ifelse( degrees, 90, pi/2 )
    r <- abs( deg90 - r )
  }
  return(r)
  }


#' Distance To A Reference Plane
#'
#' Compute the (shortest) distance between the starting point of a track and a reference plane.
#' Useful to detect directed movement and/or tracking artefacts.
#'
#' @param x a single input track; a matrix whose first column is time and whose
#'  remaining columns are a spatial coordinate.
#' @param from index, or vector of indices, of the first row of the track. If
#' \code{from} is a vector, distances are returned for all steps starting at
#' the indices in \code{from}.
#' @param xdiff row differences of x.
#' @param p1,p2,p3 numeric vectors of coordinates of three points specifying a reference plane to
#'  compute distances to.
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
distanceToPlane <- function (x, from = 1, p1 = c(0,0,0), p2 = c(0,1,0), p3 = c(1,0,0) )
{

  # Check if the given points have the correct dimensions
  plengths <- c( length(p1), length(p2), length(p3) )
  if( length(unique(plengths)) != 1){
    stop("In distanceToPlane: Points p1,p2,p3 specifying the plane must have the
         same number of coordinates.")
  }
  if( unique(plengths) != 3 ){
    stop("In distanceToPlane: Method is only supported for three-dimensional data.")
  }
  if( length(p1) != ( ncol(x) - 1 ) ){
    stop("In distanceToPlane: Plane points must have the same number of coordinates as
         coordinates in the tracking data. Currently, the method only supports 3D data.")
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

#' Distance To A Reference Point
#'
#' Compute the distance between the starting point of a track and a reference point.
#' Useful to
#' detect directed movement towards a point (see examples).
#'
#' @param x a single input track; a matrix whose first column is time and whose
#'  remaining columns are a spatial coordinate.
#' @param from index, or vector of indices, of the first row of the track. If
#' \code{from} is a vector, distances are returned for all steps starting at
#' the indices in \code{from}.
#' @param xdiff row differences of x.
#' @param p numeric vector of coordinates of the reference point p to compute distances to.
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
distanceToPoint <- function (x, from = 1, p = c(0,0,0) )
{
  # Check if the given point has the correct dimensions
  if( length(p) != ncol(x) - 1 ){
    stop("In distanceToPoint: Reference point coordinates must have the same number of dimensions as
         coordinates in the tracking data.")
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




