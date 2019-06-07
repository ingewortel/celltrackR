## ----pack, warning = FALSE, message = FALSE------------------------------
library( celltrackR )
library( ggplot2 )

## ----Tdata---------------------------------------------------------------
str( TCells, list.len = 3 )

## ------------------------------------------------------------------------
head( TCells[[1]] )

## ----bdata---------------------------------------------------------------
str( BCells, list.len = 3 )
str( Neutrophils, list.len = 3 )

## ---- fig.width = 7, fig.height=7----------------------------------------
par( mfrow = c(2,2) )
plot( TCells, main = "T cells" )
plot( BCells, main = "B cells" )
plot( Neutrophils, main = "Neutrophils" )

## ---- fig.width = 7, fig.height=7----------------------------------------
par( mfrow = c(2,2) )
plot3d( TCells, main = "T cells" )
plot3d( BCells, main = "B cells" )
plot3d( Neutrophils, main = "Neutrophils" )

## ---- fig.width=7, fig.height=7------------------------------------------
par( mfrow = c(2,2) )
plot( normalizeTracks(TCells), main = "T cells" )
plot( normalizeTracks(BCells), main = "B cells" )
plot( normalizeTracks(Neutrophils), main = "Neutrophils" )

## ------------------------------------------------------------------------
# Obtain mean speeds for each track using sapply
T.speeds <- sapply( TCells, speed )
str( T.speeds )
summary( T.speeds )

## ------------------------------------------------------------------------
hist( T.speeds )

## ------------------------------------------------------------------------
B.speeds <- sapply( BCells, speed )
N.speeds <- sapply( Neutrophils, speed )

# Create a dataframe of all data
dT <- data.frame( cells = "T cells", speed = T.speeds )
dB <- data.frame( cells = "B cells", speed = B.speeds )
dN <- data.frame( cells = "Neutrophils", speed = N.speeds )
d <- rbind( dT, dB, dN )

# Compare:
ggplot( d, aes( x = cells, y = speed ) ) +
  ggbeeswarm::geom_quasirandom() +
  scale_y_continuous( limits = c(0,NA) ) +
  theme_classic()

## ------------------------------------------------------------------------
T.steps <- subtracks( TCells, i = 1 )

## ------------------------------------------------------------------------
# Obtain instantaneous speeds for each step using sapply
T.step.speeds <- sapply( T.steps, speed )
str( T.step.speeds )
summary( T.step.speeds )

## ------------------------------------------------------------------------
hist( T.step.speeds )

## ------------------------------------------------------------------------
d.step <- data.frame( method = "step-based", speed = T.step.speeds )
d.cell <- data.frame( method = "cell-based", speed = T.speeds )
d.method <- rbind( d.step, d.cell )

ggplot( d.method, aes( x = method, y = speed ) ) +
  ggbeeswarm::geom_quasirandom( size = 0.3 ) +
  scale_y_continuous( limits = c(0,NA) ) +
  theme_classic()


## ------------------------------------------------------------------------
B.step.speeds <- sapply( subtracks(BCells,1), speed )
N.step.speeds <- sapply( subtracks(Neutrophils,1), speed )

# Create a dataframe of all data
dT <- data.frame( cells = "T cells", speed = T.step.speeds )
dB <- data.frame( cells = "B cells", speed = B.step.speeds )
dN <- data.frame( cells = "Neutrophils", speed = N.step.speeds )
dstep <- rbind( dT, dB, dN )

# Compare:
ggplot( dstep, aes( x = cells, y = speed ) ) +
  ggbeeswarm::geom_quasirandom( size = 0.3 ) +
  scale_y_continuous( limits = c(0,NA) ) +
  theme_classic()


## ------------------------------------------------------------------------
T.step.angle <- sapply( subtracks(TCells,2), overallAngle )
B.step.angle <- sapply( subtracks(BCells,2), overallAngle )
N.step.angle <- sapply( subtracks(Neutrophils,2), overallAngle )

# Create a dataframe of all data
dT <- data.frame( cells = "T cells", turn.angle = T.step.angle )
dB <- data.frame( cells = "B cells", turn.angle = B.step.angle )
dN <- data.frame( cells = "Neutrophils", turn.angle = N.step.angle )
dangle <- rbind( dT, dB, dN )

# convert radians to degrees
dangle$turn.angle <- pracma::rad2deg( dangle$turn.angle )

# Compare:
ggplot( dangle, aes( x = cells, y = turn.angle ) ) +
  ggbeeswarm::geom_quasirandom( size = 0.3 ) +
  scale_y_continuous( limits = c(0,NA) ) +
  theme_classic()


## ------------------------------------------------------------------------
x <- TCells[[8]]
x

## ------------------------------------------------------------------------
stag <- applyStaggered( x, speed, matrix = TRUE )
stag

## ------------------------------------------------------------------------
image( stag )

## ---- fig.width=4--------------------------------------------------------
filled.contour( stag )

## ---- fig.width=4--------------------------------------------------------
filled.contour( applyStaggered( TCells[[3]], speed, matrix = TRUE ) )

## ------------------------------------------------------------------------
# These are the same:
applyStaggered( TCells[[8]], speed )
mean( stag, na.rm = TRUE )

## ------------------------------------------------------------------------
# Compare the three different methods
T.stag.speed <- sapply( TCells, staggered( speed ) )
d.staggered <- data.frame( method = "staggered", speed = T.speeds )
d.method <- rbind( d.step, d.cell, d.staggered )

ggplot( d.method, aes( x = method, y = speed ) ) +
  ggbeeswarm::geom_quasirandom( size = 0.3 ) +
  scale_y_continuous( limits = c(0,NA) ) +
  theme_classic()

# cell-based and staggered cell-based are slightly different:
hist( T.speeds - T.stag.speed )

## ------------------------------------------------------------------------
aggregate( TCells, squareDisplacement, subtrack.length = 1 )

## ------------------------------------------------------------------------
aggregate( TCells, squareDisplacement, subtrack.length = 1, FUN = "mean.se" )

## ------------------------------------------------------------------------
T.msd <- aggregate( TCells, squareDisplacement, FUN = "mean.se" )
str( T.msd )

## ------------------------------------------------------------------------
ggplot( T.msd, aes( x = i, y = mean ) ) +
  geom_ribbon( aes( ymin = lower, ymax = upper) , alpha = 0.2 ,color=NA ) +
  geom_line( ) +
  labs( x = expression( paste(Delta,"t (steps)") ),
        y = "mean square displacement") +
  theme_classic()

## ------------------------------------------------------------------------
num10 <- length( subtracks( TCells, 10 ) )
num35 <- length( subtracks( TCells, 35 ) )
c( "10" = num10, "35" = num35 )

## ---- fig.width = 7------------------------------------------------------
# Combine into a single dataframe with one column indicating the celltype
# To truly compare them, report subtrack length not in number of steps but
# in their duration (which may differ between different datasets)
T.msd$cells <- "T cells"
T.msd$dt <- T.msd$i * timeStep( TCells )
B.msd <- aggregate( BCells, squareDisplacement, FUN = "mean.se" )
B.msd$cells <- "B cells"
B.msd$dt <- B.msd$i * timeStep( BCells )
N.msd <- aggregate( Neutrophils, squareDisplacement, FUN = "mean.se" )
N.msd$cells <- "Neutrophils"
N.msd$dt <- N.msd$i * timeStep( Neutrophils )
msddata <- rbind( T.msd, B.msd, N.msd )
head(msddata)

# Plot
p1 <- ggplot( msddata, aes( x = dt , y = mean, color = cells, fill = cells ) ) +
  geom_ribbon( aes( ymin = lower, ymax = upper) , alpha = 0.2 ,color=NA ) +
  geom_line( ) +
  labs( x = expression( paste(Delta,"t (seconds)") ),
        y = "mean square displacement") +
  theme_classic()
# Also make a zoomed version to look only at first part of the plot
pzoom <- p1 + coord_cartesian( xlim = c(0,500), ylim = c(0,5000) ) 
gridExtra::grid.arrange( p1, pzoom, ncol = 2 ) 

## ---- fig.width=6--------------------------------------------------------
T.acor <- aggregate( TCells, overallNormDot, FUN = "mean.se"  )
T.acor$dt <- T.acor$i * timeStep(TCells)
T.acor$cells <- "T cells"

ggplot( T.acor, aes( x = dt , y = mean, color = cells, fill = cells ) ) +
  geom_hline( yintercept = 0 ) +
  geom_ribbon( aes( ymin = lower, ymax = upper) , alpha = 0.2 ,color=NA ) +
  geom_line( ) +
  labs( x = expression( paste(Delta,"t (seconds)") ),
        y = expression(cos(alpha))) +
  theme_classic() + 
  theme( axis.line.x = element_blank() )

## ---- fig.width=6--------------------------------------------------------
# Returns TRUE if first and last step of a track are of the minimum length of 1 micron
check <- function(x){ 
  all( sapply( list(head(x,2),tail(x,2)), trackLength ) >= 1.0 )
}
# repeat analysis
T.acor2 <- aggregate( TCells, overallNormDot, 
                     FUN = "mean.se", filter.subtracks = check )
T.acor2$dt <- T.acor2$i * timeStep(TCells)
T.acor2$cells <- "T cells (filtered)"

d <- rbind( T.acor, T.acor2 )

# compare
ggplot( d, aes( x = dt , y = mean, color = cells, fill = cells ) ) +
  geom_hline( yintercept = 0 ) +
  geom_line( ) +
  labs( x = expression( paste(Delta,"t (seconds)") ),
        y = expression(cos(alpha))) +
  theme_classic() + 
  theme( axis.line.x = element_blank() )

## ---- fig.width=6--------------------------------------------------------
# autocovariance
T.acov <- aggregate( TCells, overallDot, 
                     FUN = "mean.se" )
T.acov$dt <- T.acov$i * timeStep(TCells)
T.acov$cells <- "T cells"


# compare
ggplot( T.acov, aes( x = dt , y = mean, color = cells, fill = cells ) ) +
  geom_hline( yintercept = 0 ) +
  geom_line( ) +
  labs( x = expression( paste(Delta,"t (seconds)") ),
        y = "autocovariance" ) +
  theme_classic() + 
  theme( axis.line.x = element_blank() )

## ---- fig.width = 7------------------------------------------------------
# Normalized autocovariance
T.acov[,2:4] <- T.acov[,2:4] / T.acov$mean[1]

# Other cells
B.acov <- aggregate( BCells, overallDot, FUN = "mean.se" )
B.acov$cells <- "B cells"
B.acov$dt <- B.acov$i * timeStep( BCells )
B.acov[,2:4] <- B.acov[,2:4] / B.acov$mean[1]

N.acov <- aggregate( Neutrophils, overallDot, FUN = "mean.se" )
N.acov$cells <- "Neutrophils"
N.acov$dt <- N.acov$i * timeStep( Neutrophils )
N.acov[,2:4] <- N.acov[,2:4] / N.acov$mean[1]

acovdata <- rbind( T.acov, B.acov, N.acov )

# compare
p1 <- ggplot( acovdata, aes( x = dt , y = mean, color = cells, fill = cells ) ) +
  geom_hline( yintercept = 0 ) +
  geom_ribbon( aes( ymin = lower, ymax = upper) , alpha = 0.2 ,color=NA ) +
  geom_line( ) +
  labs( x = expression( paste(Delta,"t (seconds)") ),
        y = "autocovariance" ) +
  theme_classic() + 
  theme( axis.line.x = element_blank() )

pzoom <- p1 + coord_cartesian( xlim = c(0,500) )
gridExtra::grid.arrange( p1, pzoom, ncol = 2 )

## ---- fig.width = 6, fig.height = 3--------------------------------------
par( mfrow=c(1,2) )
plot( Neutrophils )
hotellingsTest( Neutrophils, plot = TRUE, step.spacing = 3 )

## ---- fig.width = 7------------------------------------------------------
# Directional movement is 0.05 micron/sec in each dimension
directional.speed <- c( 0.05, 0.05, 0.05 )
add.dir <- function( x, speed.vector )
{
  # separate timepoints and coordinates
  tvec <- x[,1]
  coords <- x[,-1]
  
  # compute movement due to drift.
  speed.matrix <- matrix( rep( speed.vector, nrow(coords) ),
                          ncol = ncol(coords), byrow= TRUE )
  speed.matrix <- speed.matrix * tvec
  
  # Add drift to coordinates
  x[,-1] <- coords + speed.matrix
  return(x)
}

# Create data with directionality
TCells.dir <- as.tracks( lapply( TCells, add.dir, directional.speed ) )

# Plot both for comparison
par(mfrow=c(1,2) )
plot(TCells, main = "original data" )
plot(TCells.dir, main = "with directional bias" )

## ---- fig.width = 6------------------------------------------------------
step.angles <- sapply( subtracks( TCells.dir, 1), angleToDir, dvec = c(1,1,1) )
step.angles.original <- sapply( subtracks( TCells, 1 ), angleToDir, dvec = c(1,1,1) )
par(mfrow=c(1,2) )
hist( step.angles.original, main = "original data" )
hist( step.angles, main = "with directional bias" )

## ---- fig.width = 7.5, fig.height=4--------------------------------------
# Normalize the endpoint by subtracting its coordinates from all other coordinates:
normalize.endpoint <- function( x ){
  coords <- x[,-1]
  coords.norm <- t( apply( coords, 1, function(a) a - coords[nrow(coords),] ) )
  x[,-1] <- coords.norm
  return(x)
}
TCells.point <- as.tracks( lapply( TCells, normalize.endpoint ) )

# Plot for comparison:
par( mfrow = c(1,2) )
plot( TCells, main = "original data" )
plot( TCells.point, main = "with point attraction" )

## ------------------------------------------------------------------------
hotellingsTest( TCells.point, plot = TRUE, step.spacing = 3 )

## ---- fig.width = 7------------------------------------------------------
step.angles <- sapply( subtracks( TCells.point, 1), angleToPoint, p = c(1,1,1) )
step.angles.original <- sapply( subtracks( TCells, 1 ), angleToPoint, p = c(1,1,1) )
par(mfrow=c(1,2) )
hist( step.angles.original, main = "original data" )
hist( step.angles, main = "with point attraction" )

## ------------------------------------------------------------------------
step.distances <- sapply( subtracks( TCells.point, 1), distanceToPoint, p = c(1,1,1) )

df <- data.frame( dist = step.distances,
                  angles = step.angles )
ggplot( df, aes( x = dist, y = angles ) ) +
  geom_point() +
  geom_hline( yintercept = 90, color = "red" ) +
  stat_smooth( color = "black", method = "loess" ) +
  labs( x = "distance to reference point",
        y = "angle to reference point" ) +
  theme_classic()

## ------------------------------------------------------------------------
df.cells <- analyzeCellPairs( TCells )

# Plot
ggplot( df.cells, aes( x = dist, y = angle ) ) +
  geom_point( color = "gray40" ) +
  stat_smooth( span = 1, color = "black" )+
  labs( x = "distance between cell pairs",
        y = "angle between cell pairs" ) +
  geom_hline( yintercept = 90, color = "red" ) +
  theme_classic()

## ---- warning = FALSE, message = FALSE-----------------------------------
df.steps <- analyzeStepPairs( TCells, filter.steps = function(x) displacement(x)>2  )

# Plot
ggplot( df.steps, aes( x = dist, y = angle ) ) +
  geom_point( color = "gray40", size = 0.5 ) +
  stat_smooth( color = "black" )+
  labs( x = "distance between step pairs",
        y = "angle between step pairs") +
  geom_hline( yintercept = 90, color = "red" ) +
  theme_classic()


