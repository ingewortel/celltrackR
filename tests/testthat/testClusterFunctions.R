set.seed(2345)


# getFeatureMatrix
all.measures <- c(trackLength, duration, speed, displacement, squareDisplacement, maxDisplacement,
                  displacementRatio, outreachRatio, straightness, overallAngle, meanTurningAngle,
                  overallDot, overallNormDot, asphericity, hurstExponent, fractalDimension )

test_that("getFeatureMatrix returns correct output format", {
  expect_true( is.matrix( getFeatureMatrix( TCells, c( speed, meanTurningAngle ) ) ) )
  expect_equal( nrow( getFeatureMatrix( TCells, c(speed) ) ), length( TCells ) )
  expect_equal( ncol( getFeatureMatrix( TCells, c(speed) ) ), 1 )
  expect_equal( ncol( getFeatureMatrix( TCells, all.measures ) ), length(all.measures) )
  expect_equal( class( getFeatureMatrix( TCells, c( speed, meanTurningAngle ), dist = TRUE ) ), "dist" )
})
test_that("getFeatureMatrix passes arguments to dist correctly", {
  expect_equal( attr( getFeatureMatrix( TCells, c( speed, meanTurningAngle ), dist = TRUE ),"method"), "euclidean" )
  expect_equal( attr(getFeatureMatrix( TCells, c( speed, meanTurningAngle ), dist = TRUE, method = "manhattan" ),"method"), "manhattan" )
})

# clusterTracks
test_that("clusterTracks responds to input correctly", {
  expect_error( clusterTracks( TCells, c(speed), method = "hi"),
                "clusterTracks: unknown method! Please choose from: hclust, MDS, kmeans, PCA, or UMAP." )
  # kmeans requires an additional argument
  expect_error( clusterTracks( TCells, c(speed), method = "kmeans" ),
                "'centers' must be a number or a matrix" )
  expect_error( clusterTracks( TCells, c(), method = "kmeans" ),
                "clusterTracks: no measures given! Please specify at least one." )
})

test_that( "clusterTracks produces the right output", {
  # NULL returned if only plots specified, no matter the cluster method
  expect_true( is.null( clusterTracks( TCells, c(speed) ) ) )
  expect_true( is.null( clusterTracks( TCells, c(speed), method = "PCA" ) ) )
  expect_true( is.null( clusterTracks( TCells, c(speed), method = "MDS" ) ) )
  # otherwise output depends on method of choice
  expect_equal( class( clusterTracks( TCells, c(speed), method = "hclust", return.clust = TRUE ) ),
                "hclust" )
  expect_equal( class( clusterTracks( TCells, c(speed), method = "kmeans", return.clust = TRUE, centers = 2 ) ),
                "kmeans" )
  expect_equal( class( clusterTracks( TCells, c(speed), method = "PCA", return.clust = TRUE ) ),
                "matrix" )
  expect_equal( nrow( clusterTracks( TCells, c(speed), method = "PCA", return.clust = TRUE ) ),
                length(TCells) )
  expect_equal( ncol( clusterTracks( TCells, c(speed,meanTurningAngle), method = "PCA", return.clust = TRUE ) ),
                2 )
})
