# example data
im <- rjson::fromJSON( file = "immunemap.json" )
minimal.tracks <- list( list( points = list( numeric(4) ) ) )
minimal.track <- list( points = list( c(1:4), c(1:4) ) )

test_that("Importer read.immap.json checks input format correctly", {

	# Error cases:
	msg <- "Error in reading json from ImmuneMap: each track in the json should be an object that must at least contains a key 'points'. Please check json format."
	expect_error( {.read.immap.single( list() )} , msg )
	expect_error( {.read.immap.single( list( im[[1]] ) )}, msg )
	msg <- "Error in reading json from ImmuneMap: 'points' should contain an array of all numeric arrays of length 4. Your 'points' don't fit this format - please check."
	expect_error( {.read.immap.single( list( points = NULL ) )}, msg )
	expect_error( {.read.immap.single( list( points = numeric(4) ) )}, msg )
	#expect_error( {.read.immap.single( list( points = list () ))}, msg )
	msg <-  "Error in reading json from ImmuneMap: the 'points' key in the json object should contain an array of all numeric arrays of length 4. Some elements do not fulfill this criterion; please check format."
	expect_error( {.read.immap.single( list( points = list( "hi" )) )}, msg )
	expect_error( {.read.immap.single( list( points = list( numeric(3))) )}, msg )
	
	# Working case: 
	expect_is( .read.immap.single( minimal.track, warn.scaling = FALSE, keep.id = FALSE ), "tracks" )
	expect_is( read.immap.json( tracks.url="immunemap.json", scale.auto = FALSE, warn.scaling = FALSE, keep.id = FALSE, warn.celltypes = FALSE )$tracks, "tracks" )
	expect_is( read.immap.json( url="https://api.immunemap.org/video/14" )$tracks, "tracks" )
} )

test_that("Importer read.immap.json warns when scales are not set or when there are no ids to keep", {
	
	msg <- "In reading tracks from ImmuneMap: spatial scale of data unnkown, using pixels. Set parameter 'scale.pos' to supply the spatial resolution, or turn off this warning using 'warn.scaling=FALSE'."
	expect_warning( {.read.immap.single( minimal.track, scale.t = 1, keep.id = FALSE )}, msg )
	msg <-  "In reading tracks from ImmuneMap: temporal scale of data unnkown, using frames. Set parameter 'scale.t' to supply the time step between frames, or turn off this warning using 'warn.scaling=FALSE'."
	expect_warning( {.read.immap.single( minimal.track, scale.pos = 1, keep.id = FALSE )}, msg )
	msg <-  "In reading tracks from ImmuneMap json: keep.id is set to TRUE but the track contains no id. Returning a track without id. To avoid this message, set keep.id = FALSE."
	expect_warning( {.read.immap.single( minimal.track, warn.scaling = FALSE )}, msg )
	
} )

test_that("Importer read.immap.json can scale time correctly", {
	example.track <- .read.immap.single( list( points = list( c(1,2,3,4 ) ) ), warn.scaling = FALSE, keep.id = FALSE )
	example.scale.t <- 	 .read.immap.single( list( points = list( c(1,2,3,4 ) ) ), scale.t = 2, warn.scaling = FALSE, keep.id = FALSE )
	expect_true( {all( example.track[[1]][,1] * 2 == example.scale.t[[1]][,1] )} )
	for( i in 2:4 ){
		expect_true({all( example.track[[1]][,i] == example.scale.t[[1]][,i] )})
	}
} )

test_that("Importer read.immap.json can scale positions correctly", {
	example.track <- .read.immap.single( list( points = list( c(1,2,3,4 ) ) ), warn.scaling = FALSE, keep.id = FALSE )
	example.scale.pos <- 	 .read.immap.single( list( points = list( c(1,2,3,4 ) ) ), scale.pos = 2, warn.scaling = FALSE, keep.id = FALSE )
	expect_true( {all( example.track[[1]][,1] == example.scale.pos[[1]][,1] )} )
	for( i in 2:4 ){
		expect_true({all( example.track[[1]][,i] * 2 == example.scale.pos[[1]][,i] )})
	}
	scale.vec <- c(1,2,3)
	example.scale.pos2 <- .read.immap.single( list( points = list( c(1,2,3,4 ), c(1,2,3,4 ) ) ), scale.pos = scale.vec, scale.t = 1, keep.id = FALSE )
	expect_true( {all( example.track[[1]][,1] == example.scale.pos2[[1]][,1] )} )
	for( i in 2:4 ){
		expect_true({all( example.track[[1]][,i] * scale.vec[i-1] == example.scale.pos2[[1]][,i] )})
	}
} )


test_that("Importer read.immap.json returns correct output", {
	expect_is( .read.immap.single( minimal.track, keep.id = FALSE, warn.scaling = FALSE ), "tracks" )
	# Should always have x y and z coordinate, so ncol should be 4:
	expect_equal( ncol( .read.immap.single( minimal.track, keep.id = FALSE, warn.scaling = FALSE )[[1]] ), 4 )
	# check keep.id
	expect_equal( names( .read.immap.single( list( id = "hi", points = list( c(1:4), c(2:5) ) ), warn.scaling = FALSE ) ), "hi" )
} )

test_that("Function get.immap.metadata can read dates", {
	expect_false( any( is.na( read.immap.json( tracks.url="immunemap.json", scale.auto = FALSE, warn.scaling = FALSE, keep.id = FALSE, warn.celltypes = FALSE )$metadata$date ) ) )
} )
