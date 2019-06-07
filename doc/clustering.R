## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 4.5,
  fig.height = 3
)

## ----pack, warning = FALSE, message = FALSE------------------------------
library( celltrackR )

## ----Tdata---------------------------------------------------------------
str( TCells, list.len = 2 )

## ------------------------------------------------------------------------
head( TCells[[1]] )

## ----bdata---------------------------------------------------------------
str( BCells, list.len = 2 )
str( Neutrophils, list.len = 2 )

## ------------------------------------------------------------------------
T2 <- TCells
names(T2) <- paste0( "T", names(T2) )
tlab <- rep( "T", length(T2) )

B2 <- BCells
names(B2) <- paste0( "B", names(B2) )
blab <- rep( "B", length(B2) )

N2 <- Neutrophils
names(N2) <- paste0( "N", names(Neutrophils) )
nlab <- rep( "N", length( N2) )

all.tracks <- c( T2, B2, N2 )
real.celltype <- c( tlab, blab, nlab )


## ------------------------------------------------------------------------
m <- getFeatureMatrix( all.tracks, 
                       c(speed, meanTurningAngle, 
                         outreachRatio, squareDisplacement) )

# We get a matrix with a row per track and one column for each metric:
head(m)

## ---- fig.width = 4, fig.height = 3--------------------------------------
plot( m, xlab = "speed", ylab = "mean turning angle" )

## ------------------------------------------------------------------------
pca <- trackFeatureMap( all.tracks, 
               c(speed,meanTurningAngle,squareDisplacement,
                 maxDisplacement,outreachRatio ), method = "PCA", 
               labels = real.celltype, return.mapping = TRUE )

## ------------------------------------------------------------------------
pc1 <- pca[,1]
pc2 <- pca[,2]
track.speed <- sapply( all.tracks, speed )
cor.test( pc1, track.speed )

## ------------------------------------------------------------------------
trackFeatureMap( all.tracks,
               c(speed,meanTurningAngle,squareDisplacement,maxDisplacement,
                 outreachRatio ), method = "MDS",
               labels = real.celltype )

## ------------------------------------------------------------------------
trackFeatureMap( all.tracks,
        c(speed,meanTurningAngle,squareDisplacement,
          maxDisplacement,outreachRatio ), method = "UMAP",
          labels = real.celltype )

## ---- fig.width = 7, fig.height = 3--------------------------------------
clusterTracks( all.tracks,
               c(speed,meanTurningAngle,squareDisplacement,maxDisplacement,
                 outreachRatio ), method = "hclust", labels = real.celltype )

## ---- fig.height = 8, fig.width = 7--------------------------------------
clusterTracks( all.tracks,
               c(speed,meanTurningAngle,squareDisplacement,maxDisplacement,
                 outreachRatio ), method = "kmeans", 
               labels = real.celltype, centers = 3 )

