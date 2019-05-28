set.seed(2345)


# vecAngle
a <- c(1,2,3)
b <- c(1,2)
c <- matrix( c(a,a), nrow = 2, byrow = TRUE )
d <- matrix( c(a,a,a), nrow = 3, byrow = TRUE )
test_that("vecAngle responds to input correctly", {
  expect_error( vecAngle(a,c), "vecAngle: a and b must both be numeric vectors of the same length or matrices of the same size" )
  expect_error( vecAngle(a,b), "vecAngle: cannot compute angle between vectors of unequal dimensions." )
  expect_error( vecAngle(c,d), "vecAngle: cannot compute angle between vectors of unequal dimensions.")
})

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


# angleToPoint
test_that("angleToPoint responds to input correctly", {
  expect_error( angleToPoint( TCells[[1]], p = c(1,1) ),
    "In angleToPoint: Reference point coordinates must have the same number of dimensions as
         coordinates in the tracking data." )
  expect_warning( angleToPoint( TCells[[1]], p = TCells[[1]][1,-1] ),
                  "In angleToPoint: Reference point is equal to the track starting point. Angle
            is undefined, returning NA." )
})

t.steps <- subtracks( TCells, 1 )
test_that("angleToPoint returns correct output", {
  expect_true( is.numeric( angleToPoint( TCells[[1]], p = c(1,1,1) ) ) )
  expect_true( all( is.numeric( sapply( t.steps, angleToPoint, p = rnorm(3) ) ) ) )
  expect_true( all( ( sapply( t.steps, angleToPoint, p = rnorm(3) ) ) >= 0 ) )
  expect_true( all( ( sapply( t.steps, angleToPoint, p = rnorm(3) ) ) <= 180 ) )
  expect_true( all( ( sapply( t.steps, angleToPoint, p = rnorm(3), degrees = FALSE ) ) <= pi ) )
})

# angleToDir
test_that("angleToDir responds to input correctly", {
  expect_error( angleToDir( TCells[[1]], dvec = c(1,1) ),
                "In angleToDir: Direction vector must have the same number of dimensions as
         coordinates in the tracking data." )
})

test_that("angleToDir returns correct output", {
  expect_true( is.numeric( angleToDir( TCells[[1]], dvec = c(1,1,1) ) ) )
  expect_true( all( is.numeric( sapply( t.steps, angleToDir, dvec = rnorm(3) ) ) ) )
  expect_true( all( ( sapply( t.steps, angleToDir, dvec = rnorm(3) ) ) >= 0 ) )
  expect_true( all( ( sapply( t.steps, angleToDir, dvec = rnorm(3) ) ) <= 180 ) )
  expect_true( all( ( sapply( t.steps, angleToDir, dvec = rnorm(3), degrees = FALSE ) ) <= pi ) )
})

# angleToPlane
test_that("angleToPlane responds to input correctly", {
  expect_error( angleToPlane( TCells[[1]], p1 = c(1,1), p2 = c(0,0,0), p3 = c(1,1,1) ),
                "In angleToPlane: Points p1,p2,p3 specifying the plane must have the same number of coordinates." )
  expect_error( angleToPlane( projectDimensions(TCells[[1]]), p1 = c(1,1), p2 = c(0,0), p3 = c(1,0) ) ,
                "In angleToPlane: Method is only supported for three-dimensional data.")
  expect_error( angleToPlane( TCells[[1]], p1 = c(1,1), p2 = c(0,0), p3 = c(1,0) ) ,
                "In angleToPlane: Method is only supported for three-dimensional data.")
  expect_error( angleToPlane( TCells[[1]], p1 = c(1,1,1), p2 = c(0,0,0), p3 = c(1,1,1) ),
                "In angleToPlane: Points p1, p2, and p3 must be three unique points!" )
  expect_error( angleToPlane( TCells[[1]], p1 = c(1,1,1), p2 = c(0,0,0), p3 = c(2,2,2) ),
                "In angleToPlane: Points p1, p2, and p3 are on the same line and do not fully specify a plane." )
})

p1 <- rnorm(3)
p2 <- p1 + c(1,2,3)
p3 <- p1 + c(-1,2,-3)
test_that("angleToPlane returns correct output", {
  expect_true( is.numeric( angleToPlane( TCells[[1]], p1 = c(1,1,1), p2 = c(0,0,0), p3 = c(1,0,0) ) ) )
  expect_true( all( is.numeric( sapply( t.steps, angleToPlane, p1 = p1, p2 = p2, p3 = p3 ) ) ) )
  expect_true( all( ( sapply( t.steps, angleToPlane, p1 = p1, p2 = p2, p3 = p3 ) ) >= 0 ) )
  expect_true( all( ( sapply( t.steps, angleToPlane, p1 = p1, p2 = p2, p3 = p3 ) ) <= 90 ) )
  expect_true( all( ( sapply( t.steps, angleToPlane, p1 = p1, p2 = p2, p3 = p3, degrees = FALSE ) ) <= pi/2 ) )
})



# distanceToPlane
test_that("distanceToPlane responds to input correctly", {
  expect_error( distanceToPlane( TCells[[1]], p1 = c(1,1), p2 = c(0,0,0), p3 = c(1,1,1) ),
                "In distanceToPlane: Points p1,p2,p3 specifying the plane must have the same number of coordinates." )
  expect_error( distanceToPlane( TCells[[1]], p1 = c(1,1), p2 = c(0,0), p3 = c(1,0) ) ,
                "In distanceToPlane: Method is only supported for three-dimensional data.")
  expect_error( distanceToPlane( projectDimensions(TCells[[1]]), p1 = c(1,1), p2 = c(0,0), p3 = c(1,0) ) ,
                "In distanceToPlane: Method is only supported for three-dimensional data.")
  expect_error( distanceToPlane( TCells[[1]], p1 = c(1,1,1), p2 = c(0,0,0), p3 = c(1,1,1) ),
                "In distanceToPlane: Points p1, p2, and p3 must be three unique points!" )
  expect_error( distanceToPlane( TCells[[1]], p1 = c(1,1,1), p2 = c(0,0,0), p3 = c(2,2,2) ),
                "In distanceToPlane: Points p1, p2, and p3 are on the same line and do not fully specify a plane." )
})

p1 <- rnorm(3)
p2 <- p1 + c(1,2,3)
p3 <- p1 + c(-1,2,-3)
test_that("distanceToPlane returns correct output", {
  expect_true( is.numeric( distanceToPlane( TCells[[1]], p1 = c(1,1,1), p2 = c(0,0,0), p3 = c(1,0,0) ) ) )
  expect_true( all( is.numeric( sapply( t.steps, distanceToPlane, p1 = p1, p2 = p2, p3 = p3 ) ) ) )
  expect_true( all( ( sapply( t.steps, distanceToPlane, p1 = p1, p2 = p2, p3 = p3 ) ) >= 0 ) )
})

# distanceToPoint
test_that("distanceToPoint responds to input correctly", {
  expect_error( distanceToPoint( TCells[[1]], p = c(1,1) ),
                "In distanceToPoint: Reference point coordinates must have the same number of dimensions as coordinates in the tracking data." )
})

test_that("distanceToPoint returns correct output", {
  expect_true( is.numeric( distanceToPoint( TCells[[1]], p = c(1,1,1) ) ) )
  expect_true( all( is.numeric( sapply( t.steps, distanceToPoint, p = rnorm(3) ) ) ) )
  expect_true( all( ( sapply( t.steps, distanceToPoint, p = rnorm(3) ) ) >= 0 ) )
})


# angleSteps
test_that("angleSteps responds to input correctly", {
  expect_error( angleSteps( TCells, c("1","2","3"), timePoints(TCells)[1] ),
                "angleSteps: an angle is only defined for exactly 2 steps. Please provide exactly 2 trackids." )
  expect_error( angleSteps( "a", c("1","2") , timePoints(TCells)[1]),
                "angleSteps: X must be a tracks object." )
  expect_error( angleSteps( TCells[[1]], c("1","2") , timePoints(TCells)[1]),
                "angleSteps: X must be a tracks object." )
  expect_error( angleSteps( TCells, c("1","x" ), timePoints(TCells)[1] ),
                "angleSteps: cannot find all supplied trackids in the data." )
  expect_warning( angleSteps( TCells, c("1","2"), timePoints(TCells)[20] ),
                  "Warning: cannot find data for both steps. Returning NA.")
})

test_that("angleSteps returns correct output", {
  expect_true( is.numeric( angleSteps( TCells, c("1","2"), timePoints( TCells )[1] ) ) )
  expect_length( angleSteps( TCells, c("1","2"), timePoints(TCells)[1]), 1 )
  expect_true( angleSteps( TCells, c("1","2"), timePoints(TCells)[1]) >= 0 )
  expect_true( angleSteps( TCells, c("1","2"), timePoints(TCells)[1]) <= 180 )
  expect_true( angleSteps( TCells, c("1","2"), timePoints(TCells)[1], degrees = FALSE ) <= pi )
})

# distanceSteps
test_that("distanceSteps responds to input correctly", {
  expect_error( distanceSteps( TCells, c("1","2","3"), timePoints(TCells)[1] ),
                "distanceSteps: only defined for exactly 2 steps. Please provide exactly 2 trackids." )
  expect_error( distanceSteps( "a", c("1","2") , timePoints(TCells)[1]),
                "distanceSteps: X must be a tracks object." )
  expect_error( distanceSteps( TCells[[1]], c("1","2") , timePoints(TCells)[1]),
                "distanceSteps: X must be a tracks object." )
  expect_error( distanceSteps( TCells, c("1","x" ), timePoints(TCells)[1] ),
                "distanceSteps: cannot find all supplied trackids in the data." )
  expect_warning( distanceSteps( TCells, c("1","2"), timePoints(TCells)[20] ),
                  "Warning: cannot find data for both steps. Returning NA.")
})

test_that("distanceSteps returns correct output", {
  expect_true( is.numeric( distanceSteps( TCells, c("1","2"), timePoints( TCells )[1] ) ) )
  expect_length( distanceSteps( TCells, c("1","2"), timePoints(TCells)[1]), 1 )
  expect_true( distanceSteps( TCells, c("1","2"), timePoints(TCells)[1]) >= 0 )
})

# stepPairs
test_that("stepPairs responds to input correctly",{
  expect_error( stepPairs( "a"), "stepPairs: X must be a tracks object." )
  expect_error( stepPairs( TCells[[1]] ), "stepPairs: X must be a tracks object." )
})

empty.tracks <- TCells[1]
empty.tracks <- empty.tracks[-1]
single.step <- subtracks(TCells,1)[1]
test_that("stepPairs returns correct output", {
  expect_true( is.data.frame( stepPairs(TCells) ) )
  expect_true( is.data.frame( stepPairs(empty.tracks) ) )
  expect_true( is.data.frame( stepPairs(single.step) ) )
  expect_equal( ncol( stepPairs( TCells ) ), 3 )
  expect_equal( ncol( stepPairs( empty.tracks ) ), 0 )
  expect_equal( ncol( stepPairs( single.step ) ), 0 )
  expect_true( is.character( stepPairs( TCells)[,1] ) )
  expect_true( is.character( stepPairs( TCells)[,2] ) )
  expect_true( is.numeric( stepPairs( TCells )[,3] ) )
  expect_equal( ncol( stepPairs(
    TCells, filter.steps = function(t) displacement(t) < 0 ) ), 0 )
})

# analyzeStepPairs
test_that("analyzeStepPairs responds to input correctly",{
  expect_error( analyzeStepPairs( "a"), "analyzeStepPairs: X must be a tracks object." )
  expect_error( analyzeStepPairs( TCells[[1]] ), "analyzeStepPairs: X must be a tracks object." )
})

tpairs <- analyzeStepPairs( TCells )
test_that("analyzeStepPairs returns correct output", {
  expect_true( is.data.frame( tpairs ) )
  expect_true( is.data.frame( analyzeStepPairs(empty.tracks) ) )
  expect_true( is.data.frame( analyzeStepPairs(single.step) ) )
  expect_equal( ncol( tpairs ), 5 )
  expect_equal( ncol( analyzeStepPairs( empty.tracks ) ), 0 )
  expect_equal( ncol( analyzeStepPairs( single.step ) ), 0 )
  expect_true( is.character( tpairs[,1] ) )
  expect_true( is.character( tpairs[,2] ) )
  expect_true( is.numeric( tpairs[,3] ) )
  expect_true( is.numeric( tpairs[,4] ) )
  expect_true( is.numeric( tpairs[,5] ) )
  expect_equal( ncol( analyzeStepPairs(
    TCells, filter.steps = function(t) displacement(t) < 0 ) ), 0 )
})

# analyzeCellPairs
test_that("analyzeCellPairs responds to input correctly",{
  expect_error( analyzeCellPairs( "a"), "analyzeCellPairs: X must be a tracks object." )
  expect_error( analyzeCellPairs( TCells[[1]] ), "analyzeCellPairs: X must be a tracks object." )
})

tcpairs <- analyzeCellPairs( TCells )
test_that("analyzeCellPairs returns correct output", {
  expect_true( is.data.frame( tcpairs ) )
  expect_true( is.data.frame( analyzeCellPairs(empty.tracks) ) )
  expect_true( is.data.frame( analyzeCellPairs(single.step) ) )
  expect_equal( ncol( tcpairs ), 4 )
  expect_equal( ncol( analyzeCellPairs( empty.tracks ) ), 0 )
  expect_equal( ncol( analyzeCellPairs( single.step ) ), 0 )
  expect_true( is.character( tcpairs[,1] ) )
  expect_true( is.character( tcpairs[,2] ) )
  expect_true( is.numeric( tcpairs[,3] ) )
  expect_true( is.numeric( tcpairs[,4] ) )
})


# cellPairs
test_that("cellPairs responds to input correctly",{
  expect_error( cellPairs( "a"), "cellPairs: X must be a tracks object." )
  expect_error( cellPairs( TCells[[1]] ), "cellPairs: X must be a tracks object." )
})
test_that("cellPairs returns correct output", {
  expect_true( is.data.frame( cellPairs(TCells) ) )
  expect_true( is.data.frame( cellPairs(empty.tracks) ) )
  expect_true( is.data.frame( cellPairs(single.step) ) )
  expect_equal( ncol( cellPairs( TCells ) ), 2 )
  expect_equal( ncol( cellPairs( empty.tracks ) ), 0 )
  expect_equal( ncol( cellPairs( single.step ) ), 0 )
  expect_true( is.character( cellPairs( TCells)[,1] ) )
  expect_true( is.character( cellPairs( TCells)[,2] ) )
})


# angleCells
test_that("angleCells responds to input correctly",{
  expect_error( angleCells( "a", c("1","2") ), "angleCells: X must be a tracks object!" )
  expect_error( angleCells( TCells[[1]], c("1","2" ) ), "angleCells: X must be a tracks object!" )
  expect_error( angleCells( TCells, c("1","a") ), "angleCells: cannot find both cellids in data." )
})
test_that("angleCells returns correct output", {
  expect_true( is.numeric( angleCells( TCells, c("1","2" ) ) ) )
  expect_equal( angleCells( TCells, c("1","1") ), 0 )
  expect_true( angleCells( TCells, sample( names(TCells), 2 ) ) >= 0 )
  expect_true( angleCells( TCells, sample( names(TCells), 2 ) ) <= 180 )
  expect_true( angleCells( TCells, sample( names(TCells), 2 ), degrees = FALSE ) <= pi )
  expect_equal( overallAngle( subtracks( TCells, 2)[[1]], degrees = TRUE ),
                angleCells( subtracks( TCells[1], 1 ), c("0.1","0.2") ) )
})

# distanceCells
test_that("angleCells responds to input correctly",{
  expect_error( distanceCells( "a", c("1","2") ), "distanceCells: X must be a tracks object!" )
  expect_error( distanceCells( TCells[[1]], c("1","2" ) ), "distanceCells: X must be a tracks object!" )
  expect_error( distanceCells( TCells, c("1","a") ), "distanceCells: cannot find both cellids in data." )
})
dist <- distanceCells( TCells, sample( names(TCells), 2 ) )
test_that("distanceCells returns correct output", {
  expect_true( is.numeric( distanceCells( TCells, c("1","2" ) ) ) )
  expect_equal( distanceCells( TCells, c("1","1") ), 0 )
  expect_true( is.na( distanceCells( TCells, c("6","8" ) ) ) )
  expect_true(  is.na(dist) || dist >= 0 )
  expect_equal( distanceCells( subtracks( TCells, 2 ), c("0.1","0.2") ), 0 )
})
