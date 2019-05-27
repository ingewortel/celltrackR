.ulapply <- function( l, f, ... ) {
  unlist(lapply(l,f, ...), use.names=FALSE)
}

.setColAlpha <- function( col, alpha=128 ){
	do.call( grDevices::rgb, c(split(grDevices::col2rgb(col),c("red","green","blue")),
		alpha=alpha,maxColorValue=255) )
}

.minMaxOfTwoMatrices <- function(a, b, func) {
  #   print("a:")
  #   print(a)
  #   print("b:")
  #   print(b)
  res <- matrix(nrow=nrow(a), ncol=ncol(a))
  #   print(ncol(a))
  for (j in 1:ncol(a)) {
    #     print(j)
    res[1, j] <- min(a[1, j], b[1, j])
    res[2, j] <- max(a[2, j], b[2, j])
  }
  colnames(res) <- colnames(a)
  rownames(res) <- rownames(a)
  return(res)
}

.boundingBoxTrack <- function(track) {
  rbind(min = apply(track, 2, min), max = apply(track, 2, max))
}

.subtrackIndices <- function( x, overlap, lengths=nrow(x) ){
  if (overlap > x-1) {
    warning("Overlap exceeds segment length")
    overlap <- x-1
  }
  xl <- c(0,cumsum(lengths))
  r <- .ulapply(seq_along(lengths),function(l){
  	if( xl[l+1]-xl[l] <= x ){
  		c()
  	} else {
  		seq(xl[l]+1,xl[l+1]-x,x-overlap)
  	}
  })
  cbind( first=r, last=r+x )
}


.computeSegmentwiseMeans <- function(track, measure, min.segments=1, only.for.i=NULL, ...) {
  if (nrow(track) <= min.segments) {
    return(NA)
  }
  means <- c()
  debug.nr.tracks <- 0
  if (is.null(only.for.i)) {
    for (i in (min.segments):(nrow(track) - 1)) {
      subtracks <- .computeSubtracksOfISegments(track, i, ...)
      val <- sapply(subtracks, measure)
      means[i - min.segments + 1] <- mean(val)#, na.omit=TRUE)    # na.rm produces NaN if all values are NA
    }
  } else {
    if ((only.for.i >= min.segments) && (only.for.i <= (nrow(track) - 1))) {
      subtracks <- .computeSubtracksOfISegments(track, only.for.i)
      val <- sapply(subtracks, measure)
      means <- mean(val)#, na.omit=TRUE)
    } else {
      return(NA)
    }
  }
  return(means)
}

.computeSubtracksOfISegments <- function(track, i, overlap=i-1) {
  if (nrow(track) <= i) {
    return( as.tracks(list()) )
  }
  if (overlap > i-1) {
    warning("Overlap exceeds segment length")
    overlap <- i-1
  }
  l <- list()
  for (j in seq(1, (nrow(track)-i), max(1,i-overlap))) {
    l[[as.character(j)]] <- track[j:(j+i), ,drop=FALSE]
  }
  return( as.tracks(l) )
}


.segmentwiseMeans <- function(measure, ...) {
  return(function(track) {
    .computeSegmentwiseMeans(track, measure, ...)
  })
}

# Helper function of beauWalker() and is not directly called by user
# This function samples uniformly from unit sphere
.beaucheminSphereSampler <- function(d=3){
  x <- stats::rnorm(d)
  return(x/sqrt(sum(x^2)))
}

# Helper function of beauWalker() and is not directly called by user
# This function samples from the trianglular distribution
.beaucheminTriangularSampler <- function() sqrt(4*stats::runif(1))-1

# Creates a matrix that will rotate unit vector a onto unit vector b
.beaucheminRotationMatrix <- function(a,b){
		# handle undefined cases
	    theta <- acos( sum(a*b) )
    	R <- diag(rep(1,3))
		if( theta < 0.001 ){
			return(R)
		}
		if( pi-theta < 0.001 ){
			return(-R)
		}
		# compute normalized cross product of a and b
		x <- c(a[2]*b[3]-a[3]*b[2],
			a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1])
		x <- x / sqrt(sum(x^2))
		A <- rbind( c(0,-x[3],x[2]),c(x[3],0,-x[1]),c(-x[2],x[1],0) )
		R <- R + sin(theta)*A + (1-cos(theta))* (A%*%A)
    	return(R)
}

# Helper function of beauWalker() and is not directly called by user
# returns a direction for the cell to travel based on model parameters specified in beauWalker()
.beaucheminPickNewDirection <- function(
	old.direction,
	p.bias,p.persist,bias.dir,taxis.mode,
	t.free,v.free,rot.mat){
	if( stats::runif(1) < p.persist ){
		return(c(old.direction, t.free))
	}
	d <- .beaucheminSphereSampler(3)
	if( taxis.mode == 0 ){
		return(c(d,t.free))
	} else if(taxis.mode==1){
		# orthotaxis
		return(c(v.free * d*(1+p.bias*sum(d*bias.dir)),t.free))
	} else if(taxis.mode == 2){
		# topotaxis
		if( stats::runif(1) < p.bias ){
			# Approach: generate new direction as if the bias direction were (1,0,0).
			# Then rotate the resulting direction by the angles between (1,0,0) and the true
			# bias direction.
			d[1] <- .beaucheminTriangularSampler()
			circ <- .beaucheminSphereSampler(2)
			circle.scale <- sqrt(1-d[1]^2)
			d[2:3] <- circ[1:2]*circle.scale
			return(c( v.free * rot.mat%*%d, t.free))
		} else {
			return(c( v.free * d, t.free))
		}
	} else if(taxis.mode == 3){
	 	 # klinotaxis
		return(c(d*v.free,t.free*(1+p.bias*sum(d*bias.dir))))
	}
}

.gaps <- function(x, tol=0.05, deltaT=NULL){
	if( nrow(x) < 2 ){
		return(integer())
	}
	gapThreshold <- tol*deltaT
	xt <- x[,1]
	xt[apply( is.na(x), 1, any )] <-  NA
	r <- abs(diff(x[,1])-deltaT) > gapThreshold
	r[is.na(r)] <- TRUE
	which(r)
}

.normalizeTrack <- function(track) {
	cbind( track[,1,drop=FALSE], sweep(track[,-1],2,track[1,-1]) )
}

# Given a direction dir and an angle theta,
# return a random new direction (unit length)
# at angle theta to dir.
.turnByAngle2D <- function( dir, theta, degrees = TRUE ){

  if( degrees ){
    alpha <- pracma::deg2rad(alpha)
  }

  # Normalize direction of previous step
  u <- dir / sqrt( sum( dir^2) )

  # compute angle with x axis.
  alpha <- acos( u[1] )

  # New direction is theta degrees to the left or right
  if( stats::runif(1) < 0.5 ){
    phi <- alpha + theta
  } else {
    phi <- alpha - theta
  }
  v <- c( cos(phi), sin(phi) )
  v

}

# Returns a rotation matrix for 3D rotation around axis u by an angle of theta.
.rotationMatrix3D <- function( u, theta ){

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
.turnByAngle3D <- function( dir, alpha, degrees = TRUE ){

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
  ran.ang <- stats::runif(1)*2*pi
  Rm <- .rotationMatrix3D( u, ran.ang )
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

