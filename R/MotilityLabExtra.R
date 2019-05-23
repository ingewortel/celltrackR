# Additions to MotilityLab functions:
# - angle analysis
# - subsetting specific timepoints
# - finding step pairs within some distance of each other
# Some functions require dbscan (to compute nearest neighbors etc),
# or pracma (angle conversion).
require(MotilityLab)

# ===============================================
# overallAngle with degree option
# ===============================================

# Angle between two vectors a and b.
# a and b can be numeric vectors or matrices in which each row represents a numeric vector.
# In the last case, one angle is returned for each row.
vecAngle <- function( a, b, degrees = TRUE )
{
  # Checks
  if( class(a) != class(b) ){
    stop( "vecAngle: a and b must both be numeric vectors (of the same length) or matrices (of the same size)")
  }
  if( any( pracma::size(a) != pracma::size(b) ) ){
    stop( "vecAngle: cannot compute angle between vectors of unequal dimensions.")
  }
  if( is.matrix(a) && nrow(a) != nrow(b) ){
    stop( "vecAngle: a and b must have an equal number of rows.")
  }
  
  # Convert vector to matrix
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

# Adjusted overallAngle function
overallAngle <- function (x, from = 1, to = nrow(x), xdiff = diff(x), degrees = TRUE ) 
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


# ===============================================
# ANGLES/DISTANCES TO POINTS, PLANES, AND VECTORS
# ===============================================

# Compute the angle of a step with a given vector direction.
# Compute angles for all steps starting at the indices in "from".
angleToPoint <- function (x, from = 1, p = c(1,1,1), xdiff = diff(x), degrees = TRUE, quietly = FALSE )
{
  
  # Check if the given point has the correct dimensions
  if( length(p) < ncol(x) - 1 ){
    stop("In angleToPoint: Reference point coordinates must have at least as many dimensions as 
         coordinates in the tracking data.")
  } else if( (length(p) > ncol(x) - 1) & !quietly ){
    warning("In angleToPoint: Reference point coordinate set has more dimensions than coordiates
            in the tracking data. Using only the first few.")
    p <- p[1:ncol(x)]
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


# Compute the angle of a step with a given vector direction.
# Compute angles for all steps starting at the indices in "from".
angleToDir <- function (x, from = 1, dvec = c(1,1,1), xdiff = diff(x), degrees=TRUE,quietly = FALSE )
{
  
  # Check if the given direction has the correct dimensions
  if( length(dvec) < ncol(x) - 1 ){
    stop("In angleToDir: Direction vector must have at least as many dimensions as 
         coordinates in the tracking data.")
  } else if( (length(dvec) > ncol(x) - 1) & !quietly ){
    warning("In angleToDir: Direction vector has more dimensions than coordiates
            in the tracking data. Using only the first few.")
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


# Compute the angle of a step with a plane specified by three points.
# Compute angles for all steps starting at the indices in "from".
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


# Compute the distance of a step (starting point) to a plane specified by three points.
# Compute distances for all steps starting at the indices in "from"
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

# Compute the distance of a step (starting point) to a reference point.
# Compute distances for all steps starting at the indices in "from"
distanceToPoint <- function (x, from = 1, p = c(0,0,0),
                             quietly = FALSE )
{
  # Check if the given point has the correct dimensions
  if( length(p) < ncol(x) - 1 ){
    stop("In distanceToPoint: Reference point coordinates must have at least as many dimensions as 
         coordinates in the tracking data.")
  } else if( (length(p) > ncol(x) - 1) & !quietly ){
    warning("In distanceToPoint: Reference point coordinate set has more dimensions than coordiates
            in the tracking data. Using only the first few.")
    p <- p[1:ncol(x)]
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





# ===============================================
# EXTRACTING STEPS FROM GIVEN TIMEPOINT
# ===============================================

# Return all the timepoints occurring in the dataset
timePoints <- function( X )
{
  timepoints.per.track <- sapply( X, function(x) unique( x[,1] ) )
  return( unique( unlist( timepoints.per.track ) ) )
}

# Return all subtracks of i steps (i+1 positions) starting at a given timepoint t
subtracksByTime <- function( X, t, i, digits=5 )
{
  
  # Round timepoints to 4 digits; otherwise equal timepoints can seem unequal
  t <- round( t, digits )
  
  # Remove tracks that start at timepoints later than t
  # Rounding because of numerical inaccuracies (?) which make two the same
  # timepoints seem like a different number.
  X <- X[ sapply(X, function(x) round(x[1,1],digits) <= t )]
  
  # Filter out parts of the tracks at timepoints starting from t
  X2 <- as.tracks( lapply( X, function(x) x[ round(x[,1],digits) >= t, , drop=FALSE ] ) )
  
  # Get subtracks of length i starting at t by using prefixes
  sub <- prefixes( X2, i )
  
  # Remove any cells for which there are no timepoints
  sub <- sub[ !sapply(sub, is.null) ]
  
  return(as.tracks(sub))
  
}

# Subset steps from tracks with given ids at given timepoint
selectSteps <- function( X, trackids, t )
{
  X2 <- X[ as.character(trackids) ]
  return( subtracksByTime( X2, t, 1 ) )
  
}


# ===============================================
# COMPUTING DISTANCES AND ANGLES BETWEEN STEP/CELL PAIRS
# ===============================================

# Angle & distance between steps with given ids at given timepoint
angleSteps <- function( X, trackids, t, degrees = TRUE, quietly = FALSE )
{
  if( !length(trackids) == 2 ){
    stop( "angleSteps: an angle is only defined for exactly 2 steps. Please provide exactly 2 trackids.")
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
distanceSteps <- function( X, trackids, t )
{
  # Select the relevant tracks
  X2 <- X[ as.character(trackids) ]
  
  # Select the relevant timepoint, and extract coordinates (remove id/time columns)
  coords <- as.data.frame.tracks( subtracksByTime( X2, t, 0 ) )[,-c(1,2)]
  
  # Compute the euclidian distance
  diffm <- diff( as.matrix(coords) )
  return( sqrt( sum( diffm^2 ) ) )
  
}

# Return all pairs of steps in the data that:
# - occur at the same timepoint
# - fulfill a given filter criterion.
stepPairs <- function( X, filter.steps=NULL )
{
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
      pairs <- as.data.frame( t( combn( ids, 2 ) ),
                              stringsAsFactors = FALSE )
      colnames(pairs) <- c( "p1","p2" )
      
      # To a dataframe
      pairs$t <- t
      dout <- rbind( dout, pairs )
    }
  }
  dout
}

# get all pairs of cells in the data
cellPairs <- function( X )
{
  cellids <- names( X )
  pairs <- data.frame()
  
  if( length(cellids) >= 2 ){
    # Make all possible pairs of cellids
    pairs <- as.data.frame( t( combn( cellids, 2 ) ),
                            stringsAsFactors = FALSE )
    colnames(pairs) <- c( "cell1","cell2" )
  }
  return(pairs)
}

# angle between displacement vectors of two cells:
angleCells <- function( X, cellids, degrees = TRUE )
{
  X <- X[cellids]
  
  if( any( !is.element( cellids, names(X) ) ) ){
    stop( "angleCells: cannot find both cellids in data." )
  }
  
  a <- displacementVector( X[[1]] )
  b <- displacementVector( X[[2]] )
  return( vecAngle( a, b, degrees = degrees ) )
  
}

# min distance between two cells at any timepoint:
distanceCells <- function( X, cellids )
{
  
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


# Return distance and angles for ALL pairs of steps occurring at the same time (Beltman)
# Returns a dataframe with point ids, timepoint, and the angle and distance
analyzeStepPairs <- function( X, filter.steps = NULL, ... )
{
  # Obtain cell paris for each timepoint
  pairs <- stepPairs( X, filter.steps = filter.steps )
  
  # Find the distance between the step starting points
  distances <- unname( apply( pairs, 1, function(x) distanceSteps( X, x[1:2], as.numeric(x[3]) ) ) ) 
  
  # Find the angles between the steps
  angles <- unname( apply( pairs, 1, function(x) angleSteps( X, x[1:2], as.numeric(x[3]), ... ) ) ) 
  
  # Add to dataframe and return
  pairs$dist <- distances
  pairs$angle <- angles 
  pairs
}

# This function gets for each pair of tracks the angle between their
# overall displacements, as well as the minimum distance between them
# (min over all timepoints)
analyzeCellPairs <- function( X, ... )
{
  # Make all possible pairs of cellids
  pairs <- cellPairs( X )
  
  # Compute angles and distances for all cell pairs in the data
  cellangles <- apply( pairs, 1, function(x) 
    angleCells( X, x, ... ) )
  
  celldistances <- apply( pairs, 1, function(x)
    distanceCells( X, x ) )
  
  # Make a dataframe
  pairs$dist <- celldistances
  pairs$angle <- cellangles
  return(pairs)

}



# ===============================================
# CLUSTERING
# ===============================================

# Extract a feature matrix
getFeatureMatrix <- function( tracks, measures )
{
    
  values <- matrix(nrow = length(tracks))
  if (is.function(measures)) {
    measures <- c(measures)
  }
  values <- do.call(cbind, lapply(measures, 
                                  function(m) sapply(tracks, m)))

  return(values)
}

# Compute clustering
# supported methods: hclust, kmeans, cmdscale, pca, umap 
clusterTracks <- function( tracks, measures, scale = TRUE, labels = NULL, method = "hclust", return.clust = FALSE, ... ) 
{
  # Get the feature matrix
  values <- getFeatureMatrix( tracks, measures )
  
  # Scale to mean 0, sd 1 if scale = TRUE
  if (scale) {
    values <- scale(values)
  }
  
  # Compute clustering based on "method"
  if( method == "hclust" ){
    clust <- hclust( dist(values), ...)
  } else if( method == "MDS" ){
    clust <- cmdscale( dist(values), ... )
  } else if( method == "kmeans" ){
    clust <- kmeans( values, ... )
  } else if( method == "PCA" ){
    clust <- prcomp( values, ... )$x
  } else if( method == "UMAP" ){
    clust <- umap::umap( values )$layout
  } else {
    stop( "clusterTracks: unknown method! Please choose from: hclust, MDS, kmeans, PCA, or UMAP." )
  }
  
  # Make the plot
  # For dimensionality reduction methods, plot points in two dimensions,
  # colored by label if labels are given.
  if( is.element( method, c("MDS","PCA","UMAP") ) ){
      #df <- as.data.frame(clust)
      #colnames(df) <- c("V1","V2")
      if( !is.null(labels) ){
        lab <- as.numeric( factor( labels ) )
      } else {
        lab <- rep(1,length( tracks ) )
      }
      par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE )
      plot( clust, col = lab )
      if( !is.null(labels)){ 
        legend("topright", legend=unique(labels), inset=c(-0.1,0),
             col=unique(lab), pch = 1, cex=0.8)
      }
      par( mar=c(5.1, 4.1, 4.1, 2.1), xpd = FALSE )
      
  # For hierarchical clustering, plot the dendrogram colored by label if labels
      # are given.
  } else if (method == "hclust") {
    
    dend <- as.dendrogram(clust)

    # color the leaves of the dendrogram by label
    if( !is.null(labels) ){
      colors_to_use <- as.numeric( factor( labels ))
      colors_to_use <- colors_to_use[order.dendrogram(dend)]
      dendextend::labels_colors(dend) <- colors_to_use
    }
    
    plot( dend )
    if( !is.null(labels) ){
      legend("topright", legend = unique(labels),
             fill = unique( colors_to_use ), 
             border = unique( colors_to_use ), bty = "n")
    }
   
    
  # For kmeans clustering, plot each feature in the matrix for each of the k clusters.
  # Color points by label if label is given.
  } else if ( method == "kmeans" ){
    if( is.null(labels) ){
      labs <- 1
    } else {
      labs <- as.numeric( factor(labels) )
    }

    # plotting area
    nc <- 2
    nrow <- ceiling( ncol(values)/nc )
    par( mfrow=c(nrow, nc ), xpd = TRUE )
    
    # make plots for each feature
    for( i in 1:ncol(values) ){
      plot(jitter(clust$cluster), values[,i], type = "p",
                 pch = 19, col = labs,
           main = paste("Feature",i), xlab="cluster", xaxt="n")
      axis(1, at = seq(1, length(unique( clust$cluster) )))
      if( !is.null(labels) ){
        legend("topright", legend=unique(labels),
               col=unique(labs), inset=c(-0.2,0), pch = 1, cex=0.8)        
      }

    }
    
    # reset par
    par( mfrow=c(1,1), xpd = FALSE)
  }

  if( return.clust ){
    return(clust)
  }
}


# ===============================================
# SIMULATION
# ===============================================

# Given a direction dir and an angle theta,
# return a random new direction (unit length)
# at angle theta to dir.
turnByAngle2D <- function( dir, theta, degrees = TRUE ){
  
  if( degrees ){
    alpha <- pracma::deg2rad(alpha)
  }
  
  # Normalize direction of previous step
  u <- dir / sqrt( sum( dir^2) )
  
  # compute angle with x axis.
  alpha <- acos( u[1] )
  
  # New direction is theta degrees to the left or right
  if( runif(1) < 0.5 ){
    phi <- alpha + theta
  } else {
    phi <- alpha - theta
  }
  v <- c( cos(phi), sin(phi) )
  v
  
}

# Returns a rotation matrix for 3D rotation around axis u by an angle of theta.
rotationMatrix3D <- function( u, theta ){
  
  if( length(u) != 3 ){
    stop( "rotationMatrix3D: only defined for 3D coordinates." )
  }
  
  # normalize u
  u <- u / sqrt( sum( u^2) )
  
  # Build the matrix
  # See https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
  R <- matrix( 0, ncol = 3, nrow = 3 )
  for( i in 1:3 ){
    for( j in 1:3 ){
      R[i,j] <- u[i]*u[j]*(1-cos(theta))
    }
  }
  R[1,1] <- R[1,1] + cos(theta)
  R[1,2] <- R[1,2] - u[3]*sin(theta)
  R[1,3] <- R[1,3] + u[2]*sin(theta)
  R[2,1] <- R[2,1] + u[3]*sin(theta)
  R[2,2] <- R[2,2] + cos(theta)
  R[2,3] <- R[2,3] - u[1]*sin(theta)
  R[3,1] <- R[3,1] - u[2]*sin(theta)
  R[3,2] <- R[3,2] + u[1]*sin(theta)
  R[3,3] <- R[3,3] + cos(theta)
  
  return(R)
  
}

# Given a direction dir and an angle theta,
# return a random new direction (unit length)
# at angle theta to dir.
turnByAngle3D <- function( dir, alpha, degrees = TRUE ){
  
  if( degrees ){
    alpha <- pracma::deg2rad(alpha)
  }
  
  # convert direction of vector "dir" to unit vector
  # and compute angles phi and theta of spherical coordinates
  u <- dir / sqrt( sum( dir^2 ) )
  theta <- acos( u[3] )
  
  # Special case: phi is undefined if sin(theta) = 0.
  # In that case, u is directed along the z-axis so we will pick 
  # a rotation around the z-axis later. This means we can choose
  # any phi; choose zero here.
  if( sin(theta) == 0 ){
    phi <- 0
  } else {
    phi <- acos( u[1]/sin(theta) )
  }
  
  # One new vector v at angle alpha to u: keep the same phi,
  # theta2 = theta - alpha
  theta2 <- theta - alpha
  v <- c( sin(theta2)*cos(phi), 
          sin(theta2)*sin(phi),
          cos(theta2))
  
  
  # Rotate around axis of u with random angle
  ran.ang <- runif(1)*2*pi
  Rm <- rotationMatrix3D( u, ran.ang )
  vout <- as.vector( Rm %*% v )
  
  # outmatrix <- matrix(0,nrow=13, ncol=3)
  # angles <- seq(0,2*pi,length.out=12)
  # for( i in seq_along(angles) ){
  #   rm <- rotationMatrix3D( u, angles[i] )
  #   outmatrix[i+1,] <- as.vector( rm %*% v )
  # }
  # scatterplot3d( outmatrix, asp = 1, highlight.3d = TRUE )
  
  return(vout)
  
}


# Returns a simulated dataset by sampling from the speed and turning angle distributions from an 
# original track dataset.
bootstrapTrack <- function( nsteps, trackdata )
{
  
  # get dimensions of original data:
  ndim <- ncol( trackdata[[1]] ) - 1
  if( !( ndim == 2 | ndim == 3 ) ){
    stop( "bootstrapTrack: only supported in 2D or 3D.")
  }
  
  # get speeds and turning angles from original data:
  speeds <- sapply( subtracks( trackdata, 1 ), speed )
  turning.angles <- sapply( subtracks( trackdata, 2 ), overallAngle, degrees = FALSE )
  
  # get time interval dt from original data
  dt <- timeStep( trackdata )
  
  # preallocate coordinate matrix
  coords <- matrix( 0, ncol = ndim, nrow = nsteps + 1 )
  
  # loop over steps to compute coordinates
  for( s in 1:nsteps ){
    
    step.speed <- sample( speeds, 1 )
    step.disp <- dt*step.speed
    
    r <- step.disp
    
    # first step in a random direction
    if( s == 1 ){
      vec <- rnorm( ndim )
      vec.length <- sqrt( sum( vec^2 ) )
      vec <- r * vec/vec.length
      coords[s+1,] <- vec
      
    # other steps: sample turning angle from distribution
    } else {
      step.angle <- sample( turning.angles, 1 )
      
      # direction of the previous step
      u <- diff( coords[(s-1):(s),] )
      
      # compute new direction
      if( ndim == 2 ){
        new.dir <- turnByAngle2D( u, step.angle, degrees = FALSE )
      } else if ( ndim == 3 ){
        new.dir <- turnByAngle3D( u, step.angle, degrees = FALSE )
      }
      
      # multiply by displacement
      step <- new.dir * r
      
      # add to previous coord
      coords[(s+1),] <- coords[s,] + step
    }

  }
  
  m <- matrix( 0, nrow = nrow(coords), ncol = ncol(coords) + 1 )
  m[,1] <- seq(0,nrow(m)-1)*dt
  m[,-1] <- coords
  colnames(m) <- colnames( trackdata[[1]] )
  m
  
  
}
