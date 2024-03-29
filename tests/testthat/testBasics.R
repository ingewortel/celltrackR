
test_that("turning angles work",{
	expect_equal( overallAngle(
		rbind(c(0,0,0,0),c(0,1,1,1))), 0 )
	expect_equal( aggregate(TCells, overallAngle, subtrack.length=1)[1,2], 0 )
	expect_equal( aggregate(TCells, overallDot, subtrack.length=1)[1,2],
		2*mean(rowSums(do.call(rbind,normalizeTracks(subtracks(TCells,1)))[,-1]^2)) )
})

test_that("maxTrackLength works",{
	expect_equal( maxTrackLength(TCells), 40 )
	expect_equal( maxTrackLength(BCells), 40 )
	expect_equal( maxTrackLength(Neutrophils), 40 )
})

test_that("timePoints returns single vector when all tracks are equal length", {
	tp <- timePoints( TCells[c(1,1) ])
	expect_equal( ncol(tp), NULL )
} )
