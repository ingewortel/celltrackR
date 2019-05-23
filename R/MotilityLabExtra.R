# Additions to MotilityLab functions:
# - angle analysis
# - subsetting specific timepoints
# - finding step pairs within some distance of each other
# Some functions require dbscan (to compute nearest neighbors etc),
# or pracma (angle conversion).




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
