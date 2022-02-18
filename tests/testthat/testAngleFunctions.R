set.seed(2345)

## ---------- example input
# Generating inputs for the tests below.

## Datasets
TSample <- TCells[ sample( names(TCells), 10 ) ]
t.steps <- subtracks( TSample, 1 )
load( system.file("extdata", "TCellsRaw.rda", package="celltrackR" ) )
traw.sample <- TCellsRaw[ sample( names( TCellsRaw ), 10 ) ]
traw.steps <- subtracks( traw.sample, 1 )

## Plane
p1 <- rnorm(3)
p2 <- p1 + c(1,2,3)
p3 <- p1 + c(-1,2,-3)

## Specific example to test pairsByTime output. Tracks 1 and 3 do not overlap in time,
# but both do overlap with track number 2. 
test.tracks.pairsByTime <- function(){
	track.1 <- data.frame( id = "1", t = seq(1,6), x = c(1,2,2,3,3,4), y = c(1,2,3,3,2,3) )
	track.2 <- data.frame( id = "2", t = seq(3,8), x = c(3,2,3,4,4,5), y = c(5,4,4,4,5,6) )
	track.3 <- data.frame( id = "3", t = seq(7,11), x = c(5,5,6,6,7), y = c(1,2,3,4,5) )
	test.tracks <- as.tracks( rbind( track.1, track.2, track.3 ) )
}

# Make random track and a rotated counterpart such that angles are always X
rotatePoints <- function( coords, angle, degrees = TRUE ){
	if( degrees ) angle <- pracma::deg2rad(angle)
	
	single <- FALSE
	if( is.null( ncol(coords))){
		coords <- matrix( coords, ncol = 2 )
		single <- TRUE
	}
	x <- coords[,1]
	y <- coords[,2]

	
	newX <- x * cos(angle) - y * sin(angle)
	newY <- x * sin(angle) + y * cos(angle)
	
	if( single ) return( c( newX, newY ) )
	return( cbind( newX, newY ))
	
}

makeTracksAngle <- function( angle, degrees = TRUE ){
	if( degrees ) ang <- pracma::deg2rad(angle)
	track.1 <- data.frame( id = "1", t = seq( 1,5 ), x = rnorm( 5 ), y = rnorm (5) )
	track.2 <- track.1
	track.2[,c("x","y")] <- rotatePoints( track.1[,c("x","y")], angle )
	track.2[,"id"] <- "2"
	return( as.tracks( rbind( track.1, track.2 )))
}

makeTracksShifted <- function( shift ){
	track.1 <- makeRandomTrack( 5 )
	tracks <- c( track.1, track.1 )
	names(tracks) <- c("1","2") 
	tracks[["2"]][,-1] <- t( apply( tracks[[2]][,-1], 1, function(x) x + shift ) )
	return(tracks)
}

# Make random track
makeRandomTrack <- function( nCoord, id = "1", dim = 2 ){
	track <- data.frame( id = id, t = seq(1,nCoord), x = rnorm( nCoord ), y = rnorm( nCoord ) ) 
	if( dim == 3 ) track$z <- rnorm( nCoord )
	return( as.tracks( track ) )
}

# make empty tracks
makeEmptyTracks <- function(){
	return( TCells[1][-1] )
}

## ---------- vecAngle 
a <- rnorm(3)
b <- rnorm(3)
n <- 100
am <- matrix( rnorm(n*3), ncol = 3 )
bm <- matrix( rnorm(n*3), ncol = 3 )
test_that("vecAngle returns correct output type", {
  expect_true( is.numeric( vecAngle(a,b) ) )
  expect_true( is.numeric( vecAngle(am,bm) ) )
  expect_length( vecAngle( a,b ), 1 )
  expect_length( vecAngle( am,bm), n )
})
test_that("vecAngle returns degrees/radians appropriately", {
  expect_false( any( 0 > vecAngle( am, bm )) || any( vecAngle( am, bm) > 360 ) )
  expect_false( any( 0 > vecAngle( am, bm, degrees = FALSE )) || any( vecAngle( am, bm, degrees = FALSE ) > 2*pi ) )
})

a.1D <- 1
b.1D <- 4
a.2D <- c(0,1)
b.2D <- c(1,0)
a.3D <- c(0,0,1)
b.3D <- c(1,0,0)
test_that("vecAngle does not depend on input order", {
  expect_equal( vecAngle( a.1D, b.1D ), vecAngle( b.1D, a.1D ) )
  expect_equal( vecAngle( a.2D, b.2D ), vecAngle( b.2D, a.2D ) )
  expect_equal( vecAngle( a.3D, b.3D ), vecAngle( b.3D, a.3D ) )
})
test_that("vecAngle returns correct output", {
  expect_equal( vecAngle( a.1D, b.1D), 0 )
  expect_equal( vecAngle( a.1D, b.1D, degrees = FALSE), 0 )
  expect_equal( vecAngle( a.2D, b.2D ), 90 )
  expect_equal( vecAngle( a.2D, b.2D, degrees = FALSE ), pi/2 )
  expect_equal( vecAngle( a.3D, b.3D ), 90 )
  expect_equal( vecAngle( a.3D, b.3D, degrees = FALSE ), pi/2 )
})


## ---------- angleToPoint
test_that("angleToPoint returns correct output", {
  expect_true( is.numeric( angleToPoint( TCells[[1]], p = c(1,1) ) ) )
  expect_true( all( is.numeric( sapply( t.steps, angleToPoint, p = rnorm(2) ) ) ) )
  expect_true( all( ( sapply( t.steps, angleToPoint, p = rnorm(2) ) ) >= 0 ) )
  expect_true( all( ( sapply( t.steps, angleToPoint, p = rnorm(2) ) ) <= 180 ) )
  expect_true( all( ( sapply( t.steps, angleToPoint, p = rnorm(2), degrees = FALSE ) ) <= pi ) )
})

# fake track of two steps that starts at (0,0), via a random coordinate, to (0,1);
# angle to (1,0) should be 90 degrees no matter the middle coordinate.
fakeTrack <- function( D3 = FALSE ){
	fake.track <- rbind( c(0,0), 10*rnorm(2), c(0,1) )
	fake.track <- cbind( seq(1,3), fake.track )
	colnames( fake.track ) <- c("t","x","y")
	if( D3 ){
		fake.track <- cbind( fake.track, c(0, 10*rnorm(1), 0 ) )
		colnames( fake.track )[4] <- "z"
	}
	
	return( fake.track )
}
test_that("angleToPoint returns correct value", {
  expect_equal( angleToPoint( fakeTrack(), p = c(1,0) ), 90 )
})


## ---------- angleToDir
test_that("angleToDir returns correct output", {
  expect_true( is.numeric( angleToDir( TCells[[1]], dvec = c(1,1) ) ) )
  expect_true( all( is.numeric( sapply( t.steps, angleToDir, dvec = rnorm(2) ) ) ) )
  expect_true( all( ( sapply( t.steps, angleToDir, dvec = rnorm(2) ) ) >= 0 ) )
  expect_true( all( ( sapply( t.steps, angleToDir, dvec = rnorm(2) ) ) <= 180 ) )
  expect_true( all( ( sapply( t.steps, angleToDir, dvec = rnorm(2), degrees = FALSE ) ) <= pi ) )
})

# fake track of two steps that starts at (0,0), via a random coordinate, to (0,1);
# angle to direction (1,0) should be 90 degrees no matter the middle coordinate.
# function fakeTrack() defined above for angleToPoint

test_that("angleToDir returns correct value", {
  expect_equal( angleToDir( fakeTrack(), dvec = c(1,0) ), 90 )
})


## ---------- angleToPlane (plane p1,p2,p3 defined on top)
test_that("angleToPlane returns correct output", {
  expect_true( is.numeric( angleToPlane( TCellsRaw[[1]], p1 = c(1,1,1), p2 = c(0,0,0), p3 = c(1,0,0) ) ) )
  expect_true( all( is.numeric( sapply( traw.steps, angleToPlane, p1 = p1, p2 = p2, p3 = p3 ) ) ) )
  expect_true( all( ( sapply( traw.steps, angleToPlane, p1 = p1, p2 = p2, p3 = p3 ) ) >= 0 ) )
  expect_true( all( ( sapply( traw.steps, angleToPlane, p1 = p1, p2 = p2, p3 = p3 ) ) <= 90 ) )
  expect_true( all( ( sapply( traw.steps, angleToPlane, p1 = p1, p2 = p2, p3 = p3, degrees = FALSE ) ) <= pi/2 ) )
})

# fake track of two steps that starts at (0,0,0), via a random coordinate, to (0,1,0);
# angle to z-plane should be 0 degrees no matter the middle coordinate because the 
# overall displacement lies in the xy plane. Fun fakeTrack() is defined above.
test_that("angleToPlane returns correct value", {
  expect_equal( angleToPlane( fakeTrack( D3 = TRUE ), p1 = c(1,1,0), p2 = c(0,0,0), p3 = c(1,0,0) ), 0 )
})



## ---------- distanceToPlane (plane p1,p2,p3 defined on top)
test_that("distanceToPlane returns correct output", {
  expect_true( is.numeric( distanceToPlane( TCellsRaw[[1]], p1 = c(1,1,1), p2 = c(0,0,0), p3 = c(1,0,0) ) ) )
  expect_true( all( is.numeric( sapply( traw.steps, distanceToPlane, p1 = p1, p2 = p2, p3 = p3 ) ) ) )
  expect_true( all( ( sapply( traw.steps, distanceToPlane, p1 = p1, p2 = p2, p3 = p3 ) ) >= 0 ) )
})



## ---------- distanceToPoint
test_that("distanceToPoint returns correct output", {
  expect_true( is.numeric( distanceToPoint( TCells[[1]], p = c(1,1) ) ) )
  expect_true( all( is.numeric( sapply( t.steps, distanceToPoint, p = rnorm(2) ) ) ) )
  expect_true( all( ( sapply( t.steps, distanceToPoint, p = rnorm(2) ) ) >= 0 ) )
})



## ---------- pairsByTime
test_that( "pairsByTime returns correct output", {
	test.output <- pairsByTime( test.tracks.pairsByTime() )
	expect_false( any( c(1:2,9:11) %in% test.output$t ) )
	expect_true( all( 3:8 %in% test.output$t ) )
	expect_equal( test.output[ "1-2-3", "dist" ], sqrt(5) )
	expect_equal( test.output[ "1-2-4", "dist" ], sqrt(2) )
	expect_equal( test.output[ "1-2-5", "dist" ], 2 )
	expect_equal( test.output[ "1-2-6", "dist" ], 1 )
	expect_equal( test.output[ "2-3-7", "dist" ], sqrt(17) )
	expect_equal( test.output[ "2-3-8", "dist" ], 4 )
	expect_warning( output1 <- nrow( pairsByTime( TSample, searchRadius = -1 ) ),
	 "pairsByTime: no tracks share time points; returning an empty dataframe." )
	expect_equal( output1, 0 )
})



## ---------- angleSteps

test_that("angleSteps returns correct output when a single angle is requested", {
  single <- angleSteps( TCells, names(TCells)[1:2], timePoints( TCells )[2] )
  singleRadian <- angleSteps( TCells, names(TCells)[1:2], timePoints(TCells)[2], degrees = FALSE )
  expect_true( is.numeric( single ) || is.na( single ) )
  expect_length( single, 1 )
  expect_true( single >= 0 )
  expect_true( single <= 180 )
  expect_true( singleRadian <= pi )
  Angle30 <- makeTracksAngle( 30 )
  expect_equal( angleSteps( Angle30, names(Angle30)[1:2], 1 ), 30 )
})

test_that("angleSteps returns correct output when multiple angles are requested", {
  multi <- angleSteps( TCells, cbind(names(TCells)[1:3],names(TCells)[1:3]), timePoints(TCells)[2:4] )
  multiRadian <- angleSteps( TCells, cbind(names(TCells)[1:3],names(TCells)[1:3]), timePoints(TCells)[2:4], degrees = FALSE )
  expect_true( is.numeric(  multi ) )
  expect_length(  multi, 3 )
  expect_true(  all( multi >= 0 ) )
  expect_true(  all( multi <= 180 ) )
  expect_true(  all ( multiRadian <= pi ) )
  Angle30 <- makeTracksAngle( 30 )
  expect_true( all( angleSteps( Angle30, cbind(rep( names(Angle30)[1], 4),rep( names(Angle30)[2], 4)), timePoints(Angle30)[1:4] ) - 30 < 1e-10 ) )
})

test_that("angleSteps returns NA when there are no steps at the same time", {
	# tracks with ids 1 and 3 do not overlap in this track set:
	test <- test.tracks.pairsByTime()
	expect_warning( angles <- sapply( timePoints( test[1] ), function(x) angleSteps( test, c("1","3" ), x )),
		"Warning: for some pairs I cannot find data for both steps at the indicated time. Returning NA."
	)
	expect_true( all( is.na( angles )) )
})

test_that("computing angle sequentially or all at once returns the same output", {
	tr1 <- makeRandomTrack( 5 )
	tr2 <- makeRandomTrack( 5, id = "2" )
	tr <- c( tr1, tr2 )
	
	angs <- sapply( seq(1,4), function(x) angleSteps( tr, c("1","2"), x ) )
	angs2 <- angleSteps( tr, cbind( rep("1",4), rep("2",4) ), seq(1:4) )
	expect_true( all( angs == angs2 ) )
})

test_that("angles are symmetrical", {
	tr1 <- makeRandomTrack( 5 )
	tr2 <- makeRandomTrack( 5, id = "2" )
	tr <- c( tr1, tr2 )
	a1 <- angleSteps( tr, c("1","2"), 1 )
	a2 <- angleSteps( tr, c("2","1"), 1 )
	expect_equal( a1, a2 )
})

## ---------- distanceSteps
test_that("distanceSteps returns correct output when a single distance is requested", {
  single <- distanceSteps( TCells, names(TCells)[1:2], timePoints( TCells )[2] )
  expect_true( is.numeric( single ) )
  expect_length( single, 1 )
  expect_true( single >= 0 )
  shifted10x <- makeTracksShifted( c(10,0) )
  expect_equal( distanceSteps( shifted10x, names(shifted10x)[1:2], 1 ), 10 )
})

test_that("distanceSteps returns correct output when multiple distances are requested", {
  multi <- distanceSteps( TCells, cbind(names(TCells)[1:3],names(TCells)[1:3]), timePoints(TCells)[2:4] )
  multi <- distanceSteps( TCells, cbind(names(TCells)[2:4],names(TCells)[1:3]), timePoints(TCells)[2:4] )
  expect_true( is.numeric(  multi ) )
  expect_length(  multi, 3 )
  expect_true(  all( multi >= 0 ) )
  shifted10x <- makeTracksShifted( c(10,0) )
  expect_true( all( angleSteps( shifted10x, cbind(rep( names(shifted10x)[1], 4),rep( names(shifted10x)[2], 4)), timePoints(shifted10x)[1:4] ) - 30 < 1e-10 ) )
})

test_that("distanceSteps returns NA when there are no steps at the same time", {
	# tracks with ids 1 and 3 do not overlap in this track set:
	test <- test.tracks.pairsByTime()
	expect_warning( ds <- sapply( timePoints( test[1] ), function(x) distanceSteps( test, c("1","3" ), x )),
		"Warning: cannot find data for both steps. Returning NA."
	)
	expect_true( all( is.na( ds )) )
})

test_that("computing distances sequentially or all at once returns the same output", {
	tr1 <- makeRandomTrack( 5 )
	tr2 <- makeRandomTrack( 5, id = "2" )
	tr <- c( tr1, tr2 )
	
	d <- sapply( seq(1,4), function(x) distanceSteps( tr, c("1","2"), x ) )
	d2 <- distanceSteps( tr, cbind( rep("1",4), rep("2",4) ), seq(1:4) )
	expect_true( all( d == d2 ) )
})

test_that("distances are symmetrical", {
	tr1 <- makeRandomTrack( 5 )
	tr2 <- makeRandomTrack( 5, id = "2" )
	tr <- c( tr1, tr2 )
	d1 <- distanceSteps( tr, c("1","2"), 1 )
	d2 <- distanceSteps( tr, c("2","1"), 1 )
	expect_equal( d1, d2 )
})


## ---------- stepPairs

test_that("stepPairs returns correct output", {
  empty.tracks <- makeEmptyTracks()
  single.step <- subtracks( TCells[1], 1 )[1]
  p1 <- stepPairs(TSample)
  expect_true( is.data.frame( p1 ) )
  expect_true( is.data.frame( stepPairs(empty.tracks) ) )
  expect_true( is.data.frame( stepPairs(single.step) ) )
  expect_equal( ncol( p1 ), 3 )
  expect_equal( ncol( stepPairs( empty.tracks ) ), 0 )
  expect_equal( ncol( stepPairs( single.step ) ), 0 )
  expect_true( is.character( p1[,1] ) )
  expect_true( is.character( p1[,2] ) )
  expect_true( is.numeric( p1[,3] ) )
  expect_equal( ncol( stepPairs(
    TCells, filter.steps = function(t) displacement(t) < 0 ) ), 0 )
})

## ---------- analyzeStepPairs
test_that("analyzeStepPairs returns correct output", {
  tpairs <- analyzeStepPairs( TSample ) # use sample for speed
  empty.tracks <- makeEmptyTracks()
  single.step <- subtracks( TCells[1],1)[1]
  pairs.empty <- analyzeStepPairs( empty.tracks )
  pairs.single <- analyzeStepPairs( single.step )
  expect_true( is.data.frame( tpairs ) )
  expect_true( is.data.frame( pairs.empty ) )
  expect_true( is.data.frame( pairs.single ) )
  expect_equal( ncol( tpairs ), 5 )
  expect_equal( ncol( pairs.empty ), 0 )
  expect_equal( ncol( pairs.single ), 0 )
  expect_true( is.character( tpairs[,1] ) )
  expect_true( is.character( tpairs[,2] ) )
  expect_true( is.numeric( tpairs[,3] ) )
  expect_true( is.numeric( tpairs[,4] ) )
  expect_true( is.numeric( tpairs[,5] ) )
  expect_equal( ncol( analyzeStepPairs(
    TSample, filter.steps = function(t) displacement(t) < 0 ) ), 0 )
  expect_warning( noPairs <- analyzeStepPairs( TSample, searchRadius = 0 ),
    "pairsByTime: no tracks share time points; returning an empty dataframe." ) 
  expect_equal( nrow( noPairs), 0 )
})

test_that("analyzeStepPairs returns correct values", {
	test <- test.tracks.pairsByTime() 
    pbt <- pairsByTime( test )
    sp <- analyzeStepPairs( test )
    tab1 <- table( paste0( pbt[,1], "-", pbt[,2] ) )
    tab2 <- table( paste0( sp[,1], "-", sp[,2] ) )
    # For any pair in pbt, there is always one step less than there are coordinates
    expect_true( all( tab1-tab2 == 1))
	distances <- apply( sp[,1:3], 1, function(x) distanceSteps( test, c(x[1],x[2]), x[3] ) )
	angles <- apply( sp[,1:3], 1, function(x) angleSteps( test, c(x[1],x[2]), x[3] ) )
	expect_true( all( sp$dist == distances ) )
	expect_true( all( sp$angle == angles ))
})


## ---------- analyzeCellPairs

test_that("analyzeCellPairs returns correct output", {
  tcpairs <- analyzeCellPairs( TSample ) # use sample for speed
  empty.tracks <- makeEmptyTracks()
  single.step <- subtracks( TCells[1],1)[1]
  pairs.empty <- analyzeStepPairs( empty.tracks )
  pairs.single <- analyzeStepPairs( single.step )
  expect_true( is.data.frame( tcpairs ) )
  expect_true( is.data.frame( pairs.empty ) )
  expect_true( is.data.frame( pairs.single ) )
  expect_equal( ncol( tcpairs ), 4 )
  expect_equal( ncol( pairs.empty ), 0 )
  expect_equal( ncol( pairs.single ), 0 )
  expect_true( is.character( tcpairs[,1] ) )
  expect_true( is.character( tcpairs[,2] ) )
  expect_true( is.numeric( tcpairs[,3] ) )
  expect_true( is.numeric( tcpairs[,4] ) )
  # one row for each pair of ids, no matter whether they share steps or not:
  expect_equal( ncol( combn( names(TSample), 2) ), nrow( analyzeCellPairs(TSample) ) )
  expect_warning( noPairs <- analyzeCellPairs( TSample, searchRadius = 0 ),
    "pairsByTime: no tracks share time points; returning an empty dataframe." ) 
  expect_equal( nrow( noPairs), 0 )
})

test_that("analyzeCellPairs returns correct values", {
	test <- test.tracks.pairsByTime() 
    pbt <- pairsByTime( test )
    cp <- analyzeCellPairs( test )
    # in this example data, 1 and 3 do not share time points. their distance is thus
    # NA (undefined), but their angle should yield a value.
    expect_true( is.na( cp["1-3","dist" ] ) )
    expect_true( cp["1-3","angle"] >= 0 )
    expect_true( cp["1-3","angle"] <= 180 )
	distances <- apply( cp[,1:2], 1, function(x) distanceCells( test, x[1:2], quietly = TRUE ) )
	angles <- apply( cp[,1:2], 1, function(x) angleCells( test, x[1:2] ) )
	expect_true( all( cp$dist[-2] == distances[-2] ) )
	expect_true( all( cp$angle == angles ))
	expect_equal( cp$dist[1], 1 )
	expect_equal( cp$dist[3], 4 )
})


## ---------- cellPairs
test_that("cellPairs returns correct output", {
  empty.tracks <- makeEmptyTracks()
  single.step <- subtracks( TCells[1],1)[1]
  expect_true( is.data.frame( cellPairs(TCells) ) )
  expect_true( is.data.frame( cellPairs(empty.tracks) ) )
  expect_true( is.data.frame( cellPairs(single.step) ) )
  expect_equal( ncol( cellPairs( TCells ) ), 2 )
  expect_equal( ncol( cellPairs( empty.tracks ) ), 0 )
  expect_equal( ncol( cellPairs( single.step ) ), 0 )
  expect_true( is.character( cellPairs( TCells)[,1] ) )
  expect_true( is.character( cellPairs( TCells)[,2] ) )
})


## ---------- angleCells
test_that("angleCells returns correct output when a single angle is requested", {
  first <- angleCells( TCells, names(TCells)[1:2] )
  expect_true( is.numeric( first ) )
  expect_length( first, 1 )
  expect_equal( angleCells( TCells, c("1","1") ), 0 )
  sample <- angleCells( TCells, sample( names(TCells), 2 ) )
  expect_true( sample >= 0 )
  expect_true( sample <= 180 )
  expect_true( angleCells( TCells, sample( names(TCells), 2 ), degrees = FALSE ) <= pi )
  expect_equal( overallAngle( subtracks( TCells[1], 2)[[1]], degrees = TRUE ),
                angleCells( subtracks( TCells[1], 1 ), c("1.1","1.2") ) )
  Angle30 <- makeTracksAngle( 30 )
  expect_equal( angleCells( Angle30, names(Angle30)[1:2], 1 ), 30 )
})


test_that("angleCells returns correct output when multiple angles are requested", {
  multi <- angleCells( TCells, cbind(names(TCells)[2:4],names(TCells)[1:3]) )
  expect_true( is.numeric(  multi ) )
  expect_length(  multi, 3 )
  expect_true(  all( multi >= 0 ) )
  expect_true( all( multi <= 180 ) )
  Angle30 <- c( makeTracksAngle( 30 ),  makeTracksAngle( 30 ) )
  names( Angle30 ) <- c( "1","2","3","4" )
  expect_true( all( angleCells( Angle30, rbind( c("1","2"), c("3","4") ) ) - 30 < 1e-10 ) )
})

test_that("computing angles sequentially or all at once returns the same output", {
	tr <- c( makeTracksAngle( 30 ),  makeTracksAngle( 30 ) )
    names( tr ) <- c( "1","2","3","4" )
	
	a <- c( angleCells( tr, c("1","2") ) ,  angleCells( tr, c("3","4") )  )
	a2 <- angleCells( tr, rbind( c("1","2"), c("3","4") ) )
	expect_true( all( a == a2 ) )
})

test_that("angles are symmetrical", {
	tr1 <- makeRandomTrack( 5 )
	tr2 <- makeRandomTrack( 5, id = "2" )
	tr <- c( tr1, tr2 )
	a1 <- angleCells( tr, c("1","2"), 1 )
	a2 <- angleCells( tr, c("2","1"), 1 )
	expect_equal( a1, a2 )
})




## ---------- distanceCells
test_that("distanceCells returns correct output when a single distance is requested", {
  ids <- sample( names( TCells ), 2 )
  dc <- distanceCells( TCells, ids, quietly = TRUE )
  expect_true( is.numeric( dc ) || is.na(dc) )
  expect_equal( distanceCells( TCells, c("1","1") ), 0 )
  
  # NA when they share no timepts
  test <- test.tracks.pairsByTime()
  expect_true( is.na( distanceCells( test, c("1","3" ), quietly = TRUE ) ) )
  expect_true(  is.na(dc) || dc >= 0 )
  expect_equal( distanceCells( subtracks( TCells[1], 2 ), c("1.1","1.2") ), 0 )
  
  shifted10x <- makeTracksShifted( c(10,0) )
  expect_equal( distanceCells( shifted10x, names(shifted10x)[1:2] ), 10 )
  
})
test_that("distanceCells returns correct output when multiple distances are requested", {
  multi <- distanceCells( TCells, cbind(names(TCells)[2:4],names(TCells)[1:3]) )
  expect_true( is.numeric(  multi ) )
  expect_length(  multi, 3 )
  expect_true(  all( multi >= 0, na.rm = TRUE ) )
  shifted10x <- c( makeTracksShifted( c(10,0) ), makeTracksShifted( c(10,0) ) )
  names( shifted10x ) <- as.character( 1:4 )
  expect_true( all( distanceCells( shifted10x, rbind( c( "1","2"), c("3","4"))  ) - 10 < 1e-10 ) )
})

test_that("distanceCells returns NA when there are no steps at the same time", {
	# tracks with ids 1 and 3 do not overlap in this track set:
	test <- test.tracks.pairsByTime()
	expect_warning( ds <- distanceCells( test, c("1","3" )),
		"Warning: distance undefined for cells that don't share timepoints; returning NA."
	)
	expect_true( ( is.na( ds )) )
})

test_that("computing distances sequentially or all at once returns the same output", {
	tr <- c( makeTracksAngle( 30 ),  makeTracksAngle( 30 ) )
    names( tr ) <- c( "1","2","3","4" )
	
	d <- c( distanceCells( tr, c("1","2") ) ,  distanceCells( tr, c("3","4") )  )
	d2 <- distanceCells( tr, rbind( c("1","2"), c("3","4") ) )
	expect_true( all( d == d2 ) )
})

test_that("distances are symmetrical", {
	tr1 <- makeRandomTrack( 5 )
	tr2 <- makeRandomTrack( 5, id = "2" )
	tr <- c( tr1, tr2 )
	d1 <- distanceCells( tr, c("1","2") )
	d2 <- distanceCells( tr, c("2","1") )
	expect_equal( d1, d2 )
})
