load( system.file("extdata", "TCellsRaw.rda", package="celltrackR" ) )

tracks <- TCellsRaw
tracks.df <- as.data.frame.tracks( tracks )
tracks.from.df <- as.tracks.data.frame( tracks.df )

tracks2d <- projectDimensions( TCells )
tracks2d.df <- as.data.frame.tracks( tracks2d )



test_that("as.tracks.data.frame does not throw error when tracks from df with <3 dimensions", {
  expect_error(as.tracks.data.frame( tracks2d.df ), NA )
} )

test_that("as.tracks.data.frame throws an error when not at least one pos.column", {
	expect_error( as.tracks.data.frame( tracks2d.df[,1:2] ) )
} )
test_that("as.tracks.data.frame throws an error when specified pos.columns don't exist", {
	expect_error( as.tracks.data.frame( tracks2d.df, pos.columns = 3:5 ) )
} )

test_that("Tracks are converted correctly", {
  expect_equivalent(tracks[ order( names(tracks) ) ], tracks.from.df )
  expect_equivalent( as.tracks.data.frame( tracks2d.df ), tracks2d[ order( names( tracks2d ) ) ] )
} )

test_that("Tracks have correct structure", {
  expect_is(tracks.from.df[[1]], "matrix")
})

test_that("Can choose strings or factors when converting to data frame", {
	expect_s3_class(as.data.frame.tracks(TCells,idsAsFactors=TRUE)[,"id"],"factor")
	expect_type(as.data.frame.tracks(TCells,idsAsFactors=FALSE)[,"id"],"character")
})
