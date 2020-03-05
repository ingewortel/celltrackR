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
  expect_s3_class( getFeatureMatrix( TCells, c( speed, meanTurningAngle ), dist = TRUE ), "dist" )
})
test_that("getFeatureMatrix passes arguments to dist correctly", {
  expect_equal( attr( getFeatureMatrix( TCells, c( speed, meanTurningAngle ), dist = TRUE ),"method"), "euclidean" )
  expect_equal( attr(getFeatureMatrix( TCells, c( speed, meanTurningAngle ), dist = TRUE, method = "manhattan" ),"method"), "manhattan" )
})

# clusterTracks
test_that("clusterTracks responds to input correctly", {
  expect_error( clusterTracks( TCells, c(speed), method = "hi"),
                "clusterTracks: unknown method! Please choose either hclust or kmeans." )
  # kmeans requires an additional argument
  expect_error( clusterTracks( TCells, c(speed), method = "kmeans" ),
                "'centers' must be a number or a matrix" )
  expect_error( clusterTracks( TCells, c(), method = "kmeans" ),
                "clusterTracks: no measures given! Please specify at least one." )
})

test_that( "clusterTracks produces the right output", {
  # NULL returned if only plots specified, no matter the cluster method
  expect_true( is.null( clusterTracks( TCells, c(speed) ) ) )
  expect_true( is.null( clusterTracks( TCells, c(speed), method = "kmeans", centers = 3 ) ) )
  expect_true( is.null( clusterTracks( TCells, c(speed), method = "hclust" ) ) )
  # otherwise output depends on method of choice
  expect_s3_class( clusterTracks( TCells, c(speed), method = "hclust", return.clust = TRUE ),
                "hclust" )
  expect_s3_class( clusterTracks( TCells, c(speed), method = "kmeans", return.clust = TRUE, centers = 2 ),
                "kmeans" )
})


# trackFeatureMap
test_that("trackFeatureMap responds to input correctly", {
  expect_error( trackFeatureMap( TCells, c(speed), method = "hi"),
                "trackFeatureMap: unknown method! Please choose from: MDS, PCA, or UMAP." )
  expect_error( trackFeatureMap( TCells, c(), method = "PCA" ),
                "trackFeatureMap: no measures given! Please specify at least one." )
})

test_that( "trackFeatureMap produces the right output", {
  # NULL returned if only plots specified, no matter the cluster method
  expect_true( is.null( trackFeatureMap( TCells, c(speed) ) ) )
  expect_true( is.null( trackFeatureMap( TCells, c(speed), method = "MDS" ) ) )
  expect_true( is.null( trackFeatureMap( TCells, c(speed), method = "UMAP" ) ) )
  # otherwise output depends on method of choice
  expect_is( trackFeatureMap( TCells, c(speed), method = "PCA", return.mapping = TRUE ),
                "matrix" )
  expect_equal( nrow( trackFeatureMap( TCells, c(speed), method = "PCA", return.mapping = TRUE ) ),
                length(TCells) )
  # PCA returns a column (principal component) for each measure:
  expect_equal( ncol( trackFeatureMap( TCells, c(speed,meanTurningAngle), method = "PCA", return.mapping = TRUE ) ),
                2 )
  expect_equal( ncol( trackFeatureMap( TCells, c(speed), method = "PCA", return.mapping = TRUE ) ),
                1 )
  # MDS returns two columns by default but this can be tuned with k
  expect_equal( ncol( trackFeatureMap( TCells, c(speed,meanTurningAngle,straightness), method = "MDS", return.mapping = TRUE ) ),
                2 )
  expect_equal( ncol( trackFeatureMap( TCells, c(speed,meanTurningAngle,straightness), method = "MDS", return.mapping = TRUE, k = 3 ) ),
                3 )
  # UMAP returns two columns by default but this can be tuned with n_components
  expect_equal( ncol( trackFeatureMap( TCells, c(speed,meanTurningAngle,straightness), method = "UMAP", return.mapping = TRUE ) ),
                2 )
  expect_equal( ncol( trackFeatureMap( TCells, c(speed,meanTurningAngle,straightness), method = "UMAP", return.mapping = TRUE, n_components = 3 ) ),
                3 )
})

