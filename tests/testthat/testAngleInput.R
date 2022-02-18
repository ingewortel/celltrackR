set.seed(2345)
tSingle <- TCells[1]
tEmpty <- tSingle[-1]
listFormat <- lapply( TCells[1:5], function(x) x ) # this is a list but not a tracks object.

test.tracks.pairsByTime <- function(){
	track.1 <- data.frame( id = "1", t = seq(1,6), x = c(1,2,2,3,3,4), y = c(1,2,3,3,2,3) )
	track.2 <- data.frame( id = "2", t = seq(3,8), x = c(3,2,3,4,4,5), y = c(5,4,4,4,5,6) )
	track.3 <- data.frame( id = "3", t = seq(7,11), x = c(5,5,6,6,7), y = c(1,2,3,4,5) )
	test.tracks <- as.tracks( rbind( track.1, track.2, track.3 ) )
}


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


# angleToPoint
test_that("angleToPoint responds to input correctly", {
  expect_error( angleToPoint( TCells[[1]], p = c(1,1,1) ),
    "In angleToPoint: Reference point coordinates must have the same number of dimensions as
         coordinates in the tracking data." )
  expect_warning( angleToPoint( TCells[[1]], p = TCells[[1]][1,-1] ),
                  "In angleToPoint: Reference point is equal to the track starting point. Angle
            is undefined, returning NA." )
})

# angleToDir
test_that("angleToDir responds to input correctly", {
  expect_error( angleToDir( TCells[[1]], dvec = c(1,1,1) ),
                "In angleToDir: Direction vector must have the same number of dimensions as
         coordinates in the tracking data." )
})

# angleToPlane
load( system.file("extdata", "TCellsRaw.rda", package="celltrackR" ) )
test_that("angleToPlane responds to input correctly", {
  expect_error( angleToPlane( TCells[[1]], p1 = c(0,1), p2 = c(0,0,0), p3 = c(1,1) ),
                "In angleToPlane: Points p1,p2,p3 specifying the plane must have the same number of coordinates." )
  expect_error( angleToPlane( projectDimensions(TCellsRaw[[1]]), p1 = c(1,1), p2 = c(0,0), p3 = c(1,0) ) ,
                "In angleToPlane: Method is only supported for three-dimensional data.")
  expect_error( angleToPlane( TCellsRaw[[1]], p1 = c(1,1), p2 = c(0,0), p3 = c(1,0) ) ,
                "In angleToPlane: Method is only supported for three-dimensional data.")
  expect_error( angleToPlane( TCellsRaw[[1]], p1 = c(1,1,1), p2 = c(0,0,0), p3 = c(1,1,1) ),
                "In angleToPlane: Points p1, p2, and p3 must be three unique points!" )
  expect_error( angleToPlane( TCellsRaw[[1]], p1 = c(1,1,1), p2 = c(0,0,0), p3 = c(2,2,2) ),
                "In angleToPlane: Points p1, p2, and p3 are on the same line and do not fully specify a plane." )
})


# distanceToPlane
test_that("distanceToPlane responds to input correctly", {
  expect_error( distanceToPlane( TCellsRaw[[1]], p1 = c(1,1), p2 = c(0,0,0), p3 = c(1,1,1) ),
                "In distanceToPlane: Points p1,p2,p3 specifying the plane must have the same number of coordinates." )
  expect_error( distanceToPlane( TCellsRaw[[1]], p1 = c(1,1), p2 = c(0,0), p3 = c(1,0) ) ,
                "In distanceToPlane: Method is only supported for three-dimensional data.")
  expect_error( distanceToPlane( projectDimensions(TCellsRaw[[1]]), p1 = c(1,1), p2 = c(0,0), p3 = c(1,0) ) ,
                "In distanceToPlane: Method is only supported for three-dimensional data.")
  expect_error( distanceToPlane( TCellsRaw[[1]], p1 = c(1,1,1), p2 = c(0,0,0), p3 = c(1,1,1) ),
                "In distanceToPlane: Points p1, p2, and p3 must be three unique points!" )
  expect_error( distanceToPlane( TCellsRaw[[1]], p1 = c(1,1,1), p2 = c(0,0,0), p3 = c(2,2,2) ),
                "In distanceToPlane: Points p1, p2, and p3 are on the same line and do not fully specify a plane." )
})

# distanceToPoint
test_that("distanceToPoint responds to input correctly", {
  expect_error( distanceToPoint( TCells[[1]], p = c(1,1,1) ),
                "In distanceToPoint: Reference point coordinates must have the same number of dimensions as coordinates in the tracking data." )
})


# pairsByTime
test_that("pairsByTime responds to input correctly", {
	expect_error( pairsByTime( 1 ), "X must be a tracks object!" )
	expect_error( pairsByTime( "hi" ), "X must be a tracks object!" )
	expect_error( pairsByTime( TCells[[1]] ), "X must be a tracks object!" )
	expect_error( pairsByTime( listFormat ), "X must be a tracks object!" )
	expect_error( pairsByTime( TCells, times = c(48, 50) ), "times must all occur in X!" )
	expect_warning( pairsByTime( tEmpty ), "pairsByTime: tracks object X is empty. Returning an empty dataframe." )
	expect_warning( pairsByTime( tSingle ), "pairsByTime: no tracks share time points; returning an empty dataframe." )
})


# angleSteps
test_that("angleSteps responds to input correctly", {
  expect_error( angleSteps( TCells, names(TCells)[1:3], timePoints(TCells)[1] ),
                "angleSteps: an angle is only defined for exactly 2 steps. Please provide exactly 2 trackids." )
  expect_error( angleSteps( TCells, rbind(names(TCells)[1:3],names(TCells)[1:3]), timePoints(TCells)[1:2] ),
                "angleSteps: an angle is only defined for exactly 2 steps. Please provide exactly 2 trackids, or a matrix with 2 trackid columns per row." )
  expect_error( angleSteps( TCells, rbind(names(TCells)[1:2],names(TCells)[1:2]), timePoints(TCells)[1] ),
                "angleSteps: if angles are to be computed for a matrix of trackid pairs, then t must be a vector with the same length as the number of pairs." )
  expect_error( angleSteps( TCells, rbind(names(TCells)[1:2],names(TCells)[1:2]), c(0,timePoints(TCells)[1]) ),
                "angleSteps: the supplied timepoints must correspond to those in the data!" )
  expect_error( angleSteps( TCells, names(TCells)[1:2], 0 ),
                "angleSteps: the supplied timepoints must correspond to those in the data!" )
  expect_error( angleSteps( "a", names(TCells)[1:2] , timePoints(TCells)[1]),
                "angleSteps: X must be a tracks object." )
  expect_error( angleSteps( TCells[[1]], names(TCells)[1:2] , timePoints(TCells)[1]),
                "angleSteps: X must be a tracks object." )
  expect_error( angleSteps( TCells, c("1","x" ), timePoints(TCells)[1] ),
                "angleSteps: cannot find all supplied trackids in the data." )
  expect_warning( angleSteps( TCells, names(TCells)[1:2], timePoints(TCells)[20] ),
                  "Warning: for some pairs I cannot find data for both steps at the indicated time. Returning NA.")
})

# distanceSteps
test_that("distanceSteps responds to input correctly", {
  expect_error( distanceSteps( TCells, names(TCells)[1:3], timePoints(TCells)[1] ),
                "distanceSteps: only defined for exactly 2 steps. Please provide exactly 2 trackids." )
  expect_error( distanceSteps( TCells, rbind(names(TCells)[1:3],names(TCells)[1:3]), timePoints(TCells)[1:2] ),
                "distanceSteps: a distance is only defined for exactly 2 steps. Please provide exactly 2 trackids, or a matrix with 2 trackid columns per row." )
  expect_error( distanceSteps( TCells, rbind(names(TCells)[1:2],names(TCells)[1:2]), timePoints(TCells)[1] ),
                "distanceSteps: if distances are to be computed for a matrix of trackid pairs, then t must be a vector with the same length as the number of pairs." )
  expect_error( distanceSteps( TCells, rbind(names(TCells)[1:2],names(TCells)[1:2]), c(0,timePoints(TCells)[1]) ),
                "distanceSteps: the supplied timepoints must correspond to those in the data!" )
  expect_error( distanceSteps( TCells, names(TCells)[1:2], 0 ),
                "distanceSteps: the supplied timepoints must correspond to those in the data!" )

  expect_error( distanceSteps( "a", names(TCells)[1:2], timePoints(TCells)[1]),
                "distanceSteps: X must be a tracks object." )
  expect_error( distanceSteps( TCells[[1]], names(TCells)[1:2] , timePoints(TCells)[1]),
                "distanceSteps: X must be a tracks object." )
  expect_error( distanceSteps( TCells, c("1","x" ), timePoints(TCells)[1] ),
                "distanceSteps: cannot find all supplied trackids in the data." )
  expect_warning( distanceSteps( TCells, names(TCells)[1:2], timePoints(TCells)[20] ),
                  "Warning: cannot find data for both steps. Returning NA.")
})



# stepPairs
test_that("stepPairs responds to input correctly",{
  expect_error( stepPairs( "a"), "stepPairs: X must be a tracks object." )
  expect_error( stepPairs( TCells[[1]] ), "stepPairs: X must be a tracks object." )
})


# analyzeStepPairs
test_that("analyzeStepPairs responds to input correctly",{
  expect_error( analyzeStepPairs( "a"), "analyzeStepPairs: X must be a tracks object." )
  expect_error( analyzeStepPairs( TCells[[1]] ), "analyzeStepPairs: X must be a tracks object." )
})


# analyzeCellPairs
test_that("analyzeCellPairs responds to input correctly",{
  expect_error( analyzeCellPairs( "a"), "analyzeCellPairs: X must be a tracks object." )
  expect_error( analyzeCellPairs( TCells[[1]] ), "analyzeCellPairs: X must be a tracks object." )
})

# cellPairs
test_that("cellPairs responds to input correctly",{
  expect_error( cellPairs( "a"), "cellPairs: X must be a tracks object." )
  expect_error( cellPairs( TCells[[1]] ), "cellPairs: X must be a tracks object." )
})


# angleCells
test_that("angleCells responds to input correctly",{
  expect_error( angleCells( TCells, names(TCells)[1:3] ),
                "angleCells: an angle is only defined for exactly 2 cells. Please provide exactly 2 cellids." )
  expect_error( angleCells( TCells, rbind(names(TCells)[1:3],names(TCells)[1:3]), timePoints(TCells)[1:2] ),
                "angleCells: an angle is only defined for exactly 2 cells. Please provide exactly 2 cellids, or a matrix with 2 trackid columns per row." )        
  expect_error( angleCells( "a", names(TCells)[1:2]), "angleCells: X must be a tracks object!" )
  expect_error( angleCells( TCells[[1]], names(TCells)[1:2] ), "angleCells: X must be a tracks object!" )
  expect_error( angleCells( TCells, c("1","a") ), "angleCells: cannot find all cellids in data." )
  
})



# distanceCells
test_that("distanceCells responds to input correctly",{
  expect_error( distanceCells( TCells, names(TCells)[1:3] ),
                "distanceCells: a distance is only defined for exactly 2 cells. Please provide exactly 2 cellids." )
  expect_error( distanceCells( TCells, rbind(names(TCells)[1:3],names(TCells)[1:3]), timePoints(TCells)[1:2] ),
                "distanceCells: a distance is only defined for exactly 2 cells. Please provide exactly 2 cellids, or a matrix with 2 trackid columns per row." )        
  expect_error( distanceCells( "a", names(TCells)[1:2]), "distanceCells: X must be a tracks object!" )
  expect_error( distanceCells( TCells[[1]], names(TCells)[1:2] ), "distanceCells: X must be a tracks object!" )
  expect_error( distanceCells( TCells, c("1","a") ), "distanceCells: cannot find both cellids in data." )
  test <- test.tracks.pairsByTime()
  expect_warning( distanceCells( test, c("1","3") ),
                  "Warning: distance undefined for cells that don't share timepoints; returning NA.")
  expect_warning( distanceCells( test, c("1","3"), quietly = TRUE ), NA )
})
