#' Simulate an Uncorrelated Random Walk
#'
#' Generates a random track with \code{nsteps} steps in \code{dim} dimensions.
#'
#' @param nsteps desired number of steps (e.g. 10 steps generates a track with 11 positions).
#' @param dim desired number of dimensions.
#' @param mean stepwise mean drift per dimension; use 0 for an
#' unbiased Brownian motion and other values for Brownian motion with drift.
#' @param sd stepwise standard deviation per dimension.
#'
#' @details In in every step an for each dimension, a normally distributed
#' value with mean \code{mean} and standard deviation \code{sd} is
#' added to the previous cell position.
#'
#' @return A data frame  containing in cell track with \code{nsteps} steps in
#' \code{dim} dimensions is returned.
#'
#' ## The Hurst exponent of a 1D Brownian track should be near 0.5
#' hurstExponent( brownianTrack( 100, 1 ) )
#'
#' @export
brownianTrack <- function(nsteps=100, dim=3, mean=0, sd=1) {
	if( dim < 1 ){
		stop("1 or more spatial dimensions required!")
	}
	m <- cbind( 0:nsteps, replicate(dim, stats::diffinv(stats::rnorm(nsteps,mean=mean,sd=sd))))
	if( dim <= 3 ){
		colnames(m) <- c("t", "x", "y", "z")[1:(dim+1)]
	} else {
		colnames(m) <- c("t", paste0( "x", 1:dim ) )
	}
	m
}

#' Simulate a 3D Cell Track Using the Beauchemin Model
#'
#' The Beauchemin model is a simple, particle-based description of T cell motion in lymph
#' node in the absence of antigen, which is similar to a random walk (Beauchemin et al,
#' 2007).
#'
#' @param sim.time specifies the duration of the track to be generated
#' @param delta.t change in time between each timepoint.
#' @param p.persist indicates how probable a change in direction is. With p.persist = 1,
#' the direction never changes between steps and with p.persist = 0, a new direction is
#' sampled at every step.
#' @param p.bias strength of movement in the direction of \code{bias.dir}.
#' @param taxis.mode specified mode of movement. 1 := orthotaxis, 2 := topotaxis,
#' 3 := klinotaxis.
#' @param t.free time interval for how long the cell is allowed to move between turns.
#' @param v.free speed of the cell during the free motion.
#' @param t.pause time that it takes the cell to adjust movement to new direction.
#' @param bias.dir a 3D vector indicating the direction along which there is a
#' preference for movement.
#'
#' @return A track, i.e., a matrix with \code{t/delta.t} rows and 4 columns.
#'
#' @details In the Beauchemin model, cells move into a fixed direction for a fixed time \code{t.free}
#' at a fixed speed \code{v.free}. They then switch to a different direction, which is
#' sampled at uniform from a sphere. The change of direction takes a fixed time \code{t.pause},
#' during which the cell does not move. Thus, the Beauchemin model is identical to the
#' freely jointed chain model of polymer physics, except for the explicit "pause phase"
#' between subsequent steps.
#'
#' The default parameters implemented in this function
#' were found to most accurately describe 'default' T
#' cell motion in lymph nodes using least-squares fitting to the mean displacement plot
#' (Beauchemin et al, 2007).
#'
#' This function implements an extended version of the Beauchemin model, which can also
#' simulate directionally biased motion. For details, see Textor et al (2013).
#'
#' @references
#' Catherine Beauchemin, Narendra M. Dixit and Alan S. Perelson (2007), Characterizing
#' T cell movement within lymph nodes in the absence of antigen. \emph{Journal of Immunology}
#' \bold{178}(9), 5505-5512. doi:10.4049/jimmunol.178.9.5505
#'
#' Johannes Textor, Mathieu Sinn and Rob J. de Boer (2013), Analytical results on the
#' Beauchemin model of lymphocyte migration. \emph{BMC Bioinformatics} \bold{14}(Suppl 6), S10.
#' doi:10.1186/1471-2105-14-S6-S10
#'
#' @examples
#' ## Create track with model parameters and return matrix of positions
#' out <- beaucheminTrack(sim.time=20,p.persist = 0.3,taxis.mode = 1)
#' ## Plot X-Y projection
#' plot( wrapTrack(out) )
#'
#' ## Create 20 tracks and plot them all
#' out <- simulateTracks( 20, beaucheminTrack(sim.time=10,
#'   bias.dir=c(-1,1,0),p.bias=10,taxis.mode = 2,
#'   p.persist = 0.1,delta.t = 1) )
#' plot( out )
#'
#' @export
beaucheminTrack <- function(sim.time=10,delta.t=1,p.persist=0.0,p.bias=0.9,
	bias.dir=c(0,0,0),taxis.mode=1,t.free=2,v.free=18.8,t.pause=0.5){
	# Parameter checks
	if(p.persist < 0 || p.persist > 1){
		stop("p.persist must be a value between 0 and 1")
	}
	if(!(taxis.mode %in% seq(0,3))){
		stop("taxis.mode can either be 1, 2, 3, or 0 for no taxis. ",
			"See documentation for details on which model you would like to use.")
	}
	if(p.bias < 0){
		stop("p.bias must be a value greater than or equal to 0")
	}
	if(sim.time <= 0){
		stop("sim.time has to be positive")
	}
	if(delta.t <= 0){
		stop("delta.t must be a value larger than 0")
	}
	if(t.pause < 0){
		stop("t.pause must be a nonnegative value")
	}
	if(t.free <= 0){
		stop("t.free must be a value larger than 0")
	}
	if(v.free <= 0){
		stop("v.free must be a value larger than 0")
	}

	# Initializing parameters
	if(any(bias.dir != 0)){
		bias.dir <- bias.dir/sqrt(sum(bias.dir^2))
	}

	rot.mat <- NULL
	if(taxis.mode==2){
		# cache rotation matrix for topotaxis to avoid recomputing it frequently
		rot.mat <- .beaucheminRotationMatrix( c(1,0,0), bias.dir )
	}
	pos <- matrix(rep(0,4),1,4)

	d <- .beaucheminPickNewDirection( NULL,p.bias,0,bias.dir,taxis.mode,
			t.free,v.free,rot.mat )
	p <- c(0,0,0)
	t <- t.pause
	while( t <= sim.time+t.free+t.pause ){
		pnew <- p+d[4]*d[-4]
		tnew <- t+d[4]
		pos <- rbind( pos, c(t,p), c(tnew,pnew) )
		t <- tnew + t.pause
		p <- pnew
		if( p.persist == 0 || stats::runif(1) > p.persist ){
			d <- .beaucheminPickNewDirection( d,p.bias,0,bias.dir,taxis.mode,
					t.free,v.free,rot.mat )
		}
	}
	pnew <- p+d[4]*d[-4]
	tnew <- t+d[4]
	pos <- rbind( pos, c(t,p), c(tnew,pnew) )


	# interpolate track observations according to delta.t
	t <- seq(0,sim.time,by=delta.t)+stats::runif(1,max=t.free+t.pause)
	pos.interpolated <- apply(pos[,-1],2,function(x) stats::approx(pos[,1], x, xout=t)$y)
	pos <- cbind( t, pos.interpolated )
	colnames(pos) <- c("t","x","y","z")

	return(.normalizeTrack(pos))
}


#' Generate Tracks by Simulation
#'
#' Generic function that executes \code{expr}, which is expected to
#' return a track, \code{n} times and stores the output in a \code{tracks}
#' object. Basically, this works like \code{\link{replicate}} but for tracks.
#'
#' @param n number of tracks to be generated.
#' @param expr the expression, usually a call, that generates a single track.
#'
#' @return A \code{tracks} object containing \code{n} tracks.
#'
#' @examples
#' ## Generate 10 tracks, 100 steps each, from a random walk with standard normally
#' ## distributed increments and plot them
#' plot( simulateTracks( 10, brownianTrack(100,3) ) )
#'
#' @export
simulateTracks <- function( n, expr ){
	as.tracks(
		sapply( as.character(seq_len(n)), eval.parent(substitute(function(...) expr)),
		simplify=FALSE, USE.NAMES=TRUE ) )
}

#' Simulate Tracks via Bootstrapping of Speed and Turning Angle from a Real Track Dataset
#'
#' Returns a simulated dataset by sampling from the speed and turning angle distributions from an
#' original track dataset (only in 2 or 3 dimensions)
#'
#' @param nsteps desired number of steps (e.g. 10 steps generates a track with 11 positions).
#' @param trackdata a tracks object to extract speeds and turning angles from.
#'
#' @details The number of dimensions is kept the same as in the original data (if data is 3D but
#' simulated tracks should be 2D, consider calling \code{\link{projectDimensions}} on the input
#' data before supplying it to \code{bootstrapTrack}). The time interval between "measurements"
#' of the simulated track equals that in the real data and is found via \code{\link{timeStep}}.
#' The first step starts at the origin in a random direction, with a speed sampled from the
#' speed distribution to determine its displacement. All subsequent steps also have their
#' turning angles sampled from the turning angle distribution in the data.
#'
#' @return A data frame  containing in cell track with \code{nsteps} steps in
#' the same number of dimensions as the original data is returned.
#'
#' ## Generate bootstrapped tracks of the TCell data; compare its speed distribution to the
#' ## original data (should be the same).
#' T.bootstrap <- bootstrapTrack( 100, TCells )
#' step.speeds.real <- sapply( subtracks(TCells,1), speed )
#' step.speeds.bootstrap <- sapply( subtracks( T.bootstrap, 1), speed )
#' qqplot( step.speeds.real, step.speeds.bootstrap )
#'
#'
#' @export
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
      vec <- stats::rnorm( ndim )
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
        new.dir <- .turnByAngle2D( u, step.angle, degrees = FALSE )
      } else if ( ndim == 3 ){
        new.dir <- .turnByAngle3D( u, step.angle, degrees = FALSE )
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
