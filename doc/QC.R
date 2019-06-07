## ----pack, warning = FALSE, message = FALSE------------------------------
library( celltrackR )
library( ggplot2 )

## ----Tdata---------------------------------------------------------------
str( TCells, list.len = 3 )

## ------------------------------------------------------------------------
head( TCells[[1]] )

## ----tracklength-check---------------------------------------------------
# Each track has a coordinate matrix with one row per coordinate;
# The number of steps is the number of rows minus one.
track.lengths <- sapply( TCells, nrow ) - 1
hist( track.lengths, xlab = "Track length (#steps)" )
summary( track.lengths )

## ----max-tracklength-----------------------------------------------------
# This is the number of coordinates, so the number of steps is one less.
maxTrackLength( TCells )

## ----filter-short--------------------------------------------------------
# nrow() of a track is always the number of steps plus one.
# For steps >= min.steps, we can substitute nrow > min.steps:
filter.criterion <- function( x, min.steps ){
  nrow(x) > min.steps
}

TCells.filtered <- filterTracks( filter.criterion, TCells, min.steps = 11 )
# Or shorthand: filterTracks( function(x) nrow(x) > min.steps, TCells )

## ----check-filtered, fig.width=6-----------------------------------------
# Find lengths of the filtered dataset
track.lengths.filtered <- sapply( TCells.filtered, nrow ) - 1

# Histograms of track lengths before and after
par( mfrow=c(1,2) )
hist( track.lengths, xlab = "Track length (#steps)",
      main = "Before filtering", breaks = seq(0, 40, by = 10 ) )
hist( track.lengths.filtered, xlab = "Track length (#steps)",
      main = "After filtering",  breaks = seq(0, 40, by = 10 ))

# Check how many tracks are left:
length( TCells.filtered )

## ---- fig.width=7--------------------------------------------------------
# Drift is 0.05 micron/sec in each dimension
drift.speed <- c( 0.05, 0.05, 0.05 )
add.drift <- function( x, drift.vector )
{
  # separate timepoints and coordinates
  tvec <- x[,1]
  coords <- x[,-1]
  
  # compute movement due to drift.
  drift.matrix <- matrix( rep( drift.vector, nrow(coords) ),
                          ncol = ncol(coords), byrow= TRUE )
  drift.matrix <- drift.matrix * tvec
  
  # Add drift to coordinates
  x[,-1] <- coords + drift.matrix
  return(x)
}

# Create data with drift
TCells.drift <- as.tracks( lapply( TCells, add.drift, drift.speed ) )

# Plot both for comparison
par(mfrow=c(1,2) )
plot(TCells, main = "original data" )
plot(TCells.drift, main = "with tissue drift" )

# Overlay starting points to view directionality
plot( normalizeTracks(TCells), main = "original data" )
plot( normalizeTracks(TCells.drift), main = "with tissue drift" )


## ------------------------------------------------------------------------
hotellingsTest( TCells, plot = TRUE, col = "gray" )

## ------------------------------------------------------------------------
hotellingsTest( TCells, plot = TRUE, col = "gray", step.spacing = 5 )

## ------------------------------------------------------------------------
hotellingsTest( TCells.drift, plot = TRUE, col = "gray", step.spacing = 5 )

## ------------------------------------------------------------------------
hotellingsTest( TCells, dim = c("x","y","z"), step.spacing = 5 )
hotellingsTest( TCells.drift, dim = c("x","y","z"), step.spacing = 5 )

## ---- echo = FALSE, fig.width=7------------------------------------------
# sp <- seq(0,10)
# 
# htest.nodrift <- lapply( sp, function(x) hotellingsTest( TCells, 
#                                                  dim = c("x","y","z"),
#                                                  step.spacing = x ))
# htest.drift <- lapply( sp, function(x) hotellingsTest( TCells.drift, 
#                                                  dim = c("x","y","z"),
#                                                  step.spacing = x ))
# 
# pval.nodrift <- sapply( htest.nodrift, function(x) x$p.value )
# pval.drift <- sapply( htest.drift, function(x) x$p.value )
# stat.nodrift <- sapply( htest.nodrift, function(x) unname(x$statistic) )
# stat.drift <- sapply( htest.drift, function(x) unname(x$statistic) )
# 
# d.drift <- data.frame( step.spacing = sp,
#                        exp = "drift",
#                        pval = pval.drift,
#                        T.statistic = stat.drift )
# d.nodrift <- data.frame( step.spacing = sp,
#                        exp = "no drift",
#                        pval = pval.nodrift,
#                        T.statistic = stat.nodrift )
# d <- rbind( d.drift, d.nodrift )
# 
# p1 <- ggplot( d, aes( x = step.spacing, y = pval, color = exp ) ) +
#   geom_line() +
#   geom_hline( yintercept = 0.05, color = "red" ) +
#   scale_color_manual( values = c( "drift" = "black", "no drift" = "gray")  ) +
#   scale_y_log10() +
#   theme_classic()
# 
# p2 <- ggplot( d, aes( x = step.spacing, y = T.statistic, color = exp ) ) +
#   geom_line() +
#   scale_y_continuous( limits = c(0,NA) ) +
#   scale_color_manual( values = c( "drift" = "black", "no drift" = "gray")  ) +
#   theme_classic()
# 
# gridExtra::grid.arrange(p1,p2,ncol=2)

## ---- warning = FALSE, message = FALSE, fig.width=7----------------------
# compute for both original as drift data
df.drift <- analyzeCellPairs( TCells.drift )
df.norm <- analyzeCellPairs( TCells )

# Plot
p.norm <- ggplot( df.norm, aes( x = dist, y = angle ) ) +
  geom_point( color = "gray40" ) +
  stat_smooth( span = 1, color = "black" )+
  labs( x = "distance between cell pairs",
        y = "angle between cell pairs",
        title = "original data") +
  geom_hline( yintercept = 90, color = "red" ) +
  theme_classic()

p.drift <- ggplot( df.drift, aes( x = dist, y = angle ) ) +
  geom_point( color = "gray40" ) +
  stat_smooth( span = 1, color = "black" )+
  labs( x = "distance between cell pairs",
        y = "angle between cell pairs",
        title = "data with drift") +
  geom_hline( yintercept = 90, color = "red" ) +
  theme_classic()

gridExtra::grid.arrange( p.norm, p.drift, ncol = 2 )

## ---- warning = FALSE, message = FALSE, fig.width=7----------------------
# compute for both original as drift data
df.drift <- analyzeStepPairs( TCells.drift, filter.steps = function(x) displacement(x)>2  )
df.norm <- analyzeStepPairs( TCells, filter.steps = function(x) displacement(x)>2  )

# Plot
p.norm <- ggplot( df.norm, aes( x = dist, y = angle ) ) +
  geom_point( color = "gray40", size = 0.5 ) +
  stat_smooth( span = 0.75, color = "black" )+
  labs( x = "distance between step pairs",
        y = "angle between step pairs",
        title = "original data") +
  geom_hline( yintercept = 90, color = "red" ) +
  theme_classic()

p.drift <- ggplot( df.drift, aes( x = dist, y = angle ) ) +
  geom_point( color = "gray40", size = 0.5 ) +
  stat_smooth( span = 0.75, color = "black" )+
  labs( x = "distance between step pairs",
        y = "angle between step pairs",
        title = "data with drift") +
  geom_hline( yintercept = 90, color = "red" ) +
  theme_classic()

gridExtra::grid.arrange( p.norm, p.drift, ncol = 2 )

## ------------------------------------------------------------------------
# Get steps and find their displacement vectors
steps.drift <- subtracks( TCells.drift, 1 )
step.disp <- t( sapply( steps.drift, displacementVector ) )

# Get the mean
mean.displacement <- colMeans( step.disp )

# Divide this by the mean timestep to get a drift speed
drift.speed <- mean.displacement/timeStep( TCells.drift )
drift.speed

## ---- fig.width = 7------------------------------------------------------
correct.drift <- function( x, drift.vector )
{
  # separate timepoints and coordinates
  tvec <- x[,1]
  coords <- x[,-1]
  
  # compute movement due to drift.
  drift.matrix <- matrix( rep( drift.vector, nrow(coords) ),
                          ncol = ncol(coords), byrow= TRUE )
  drift.matrix <- drift.matrix * tvec
  
  # Add drift to coordinates
  x[,-1] <- coords - drift.matrix
  return(x)
}

# Create data with drift
TCells.corrected <- as.tracks( lapply( TCells.drift, correct.drift, drift.speed ) )

# Compare
par( mfrow = c(1,2) )
plot( TCells, col = "gray", main = "uncorrected" )
plot( TCells.drift, col = "red", add = TRUE )
plot( TCells, col = "gray", main = "corrected" )
plot( TCells.corrected, col = "red", add = TRUE )

#
par( mfrow = c(1,2) )
plot( normalizeTracks(TCells), col = "gray", main = "uncorrected" )
plot( normalizeTracks(TCells.drift), col = "red", add = TRUE )
plot( normalizeTracks(TCells), col = "gray", main = "corrected" )
plot( normalizeTracks(TCells.corrected), col = "red", add = TRUE )

## ------------------------------------------------------------------------
# Take the track with id "2"
dup.track <- TCells[["2"]]

# Add some noise to coordinates
dup.track[,"x"] <- dup.track[,"x"] + rnorm( nrow(dup.track), sd = 0.5 )
dup.track[,"y"] <- dup.track[,"y"] + rnorm( nrow(dup.track), sd = 0.5 )
dup.track[,"z"] <- dup.track[,"z"] + rnorm( nrow(dup.track), sd = 0.5 )

# Wrap the track in a tracks object and add it to the TCell data with
# a unique id number
dup.track <- wrapTrack( dup.track )
names(dup.track) <- "22"
TCells.dup <- c( TCells, dup.track )

## ---- warning = FALSE----------------------------------------------------
df <- analyzeCellPairs( TCells.dup )

# label cellpairs that have both angle and distance below threshold
angle.thresh <- 90 # in degrees
dist.thresh <- 10 # this should be the expected cell radius
df$id <- paste0( df$cell1,"-",df$cell2 )
df$id[ !(df$angle < angle.thresh & df$dist < dist.thresh) ] <- "" 

# Plot
ggplot( df, aes( x = dist, y = angle ) ) +
  geom_point( color = "gray40" ) +
  geom_text( aes( label = id ), color = "red" ) +
  labs( x = "distance between cell pairs",
        y = "angle between cell pairs" ) +
  geom_hline( yintercept = angle.thresh, col = "blue",lty=2 ) +
  geom_vline( xintercept = dist.thresh, col = "blue", lty=2) +
  theme_classic()


## ------------------------------------------------------------------------
plot( TCells.dup[c("2","22")])

## ------------------------------------------------------------------------
tracks <- TCells
bb <- boundingBox( tracks )
bb

## ------------------------------------------------------------------------
# Define points:
lower1 <- c( bb["min","x"], bb["min","y"], bb["min","z"] )
lower2 <- c( bb["max","x"], bb["min","y"], bb["min","z"] )
lower3 <- c( bb["max","x"], bb["max","y"], bb["min","z"] )
zsize <- bb["max","z"] - bb["min","z"]

# Compute angles and distances of steps to this plane.
single.steps <- subtracks( tracks, 1 )
angles <- sapply( single.steps, angleToPlane, p1 = lower1, 
                  p2 = lower2, p3 = lower3 )
distances <- sapply( single.steps, distanceToPlane, p1 = lower1,
                     p2 = lower2, p3 = lower3 )
df <- data.frame( angles = angles,
                  distances = distances )

# Plot
ggplot( df, aes( x = distances,
                 y = angles ) ) +
  geom_point( color = "gray40" ) +
  stat_smooth( method = "loess", span = 1, color = "black" ) +
  geom_hline( yintercept = 32.7, color = "red" ) +
  scale_x_continuous( limits=c(0,zsize ), expand = c(0,0) ) +
  theme_classic()

## ----get-steps-----------------------------------------------------------
# Extract all subtracks of length 1 (that is, all "steps")
single.steps <- subtracks( TCells, 1 )

# The output is a new tracks object with a unique track for each step
# in the data (no longer grouped by the original cell they came from):
str( single.steps, list.len = 3 )

## ----check-avdt----------------------------------------------------------
median.dt <- timeStep( TCells )
median.dt

## ----check-dt------------------------------------------------------------
step.dt <- sapply( single.steps, duration )
str(step.dt)

## ----check-dt-hist-------------------------------------------------------
dt.diff.perc <- (step.dt - median.dt) * 100 / median.dt
hist( dt.diff.perc, xlab = "dt (percentage difference from median)" )

## ------------------------------------------------------------------------
# This function randomly removes coordinates from a track dataset with probability "prob"
remove.points <- function( track, prob=0.1 ){
  
    tlength <- nrow( track )
    remove.rows <- sample( c(TRUE,FALSE), tlength, replace=TRUE,
                           prob = c(prob, (1-prob) ) )
    track <- track[!remove.rows,]
  return(track)
}

# Apply function to dataset to randomly remove coordinates in the data
TCells.gap <- as.tracks( lapply( TCells, remove.points ) )

## ------------------------------------------------------------------------
# median dt of the new data
median.dt.gap <- timeStep( TCells.gap )

# duration of the individual steps
steps.gap <- subtracks( TCells.gap, 1 )
step.dt.gap <- sapply( steps.gap, duration )

# express difference as percentage of median dt
dt.diff.perc.gap <- (step.dt.gap - median.dt.gap) * 100 / median.dt.gap
hist( dt.diff.perc.gap, xlab = "dt (percentage difference from median)" )

## ------------------------------------------------------------------------

T.step.disp <- sapply( single.steps, displacement )
T.gap.disp <- sapply( steps.gap, displacement )

lapply( list( original = T.step.disp, gaps = T.gap.disp ), summary )

## ------------------------------------------------------------------------
T.norm.disp <- sapply( single.steps, normalizeToDuration( displacement ) )
T.norm.gap.disp <- sapply( steps.gap, normalizeToDuration( displacement ) )
lapply( list( original = T.norm.disp, gaps = T.norm.gap.disp ), summary )

## ---- fig.width=6--------------------------------------------------------

# Repair gaps by splitting or interpolation, the number of tracks is 
# different after each fix
split.gap <- repairGaps( TCells.gap, how = "split" )
interpolate.gap <- repairGaps( TCells.gap, how = "interpolate" )

c( "after splitting" = length( split.gap),
   "after interpolation" = length( interpolate.gap ) )

## ------------------------------------------------------------------------
T2 <- subsample( TCells, k = 2 )

## ------------------------------------------------------------------------
# displacement
T.steps <- subtracks( TCells, 1 )
T.disp <- sapply( T.steps, displacement )

T2.steps <- subtracks( T2, 1 )
T2.disp <- sapply( T2.steps, displacement )

lapply( list( T = T.disp, T2 = T2.disp ), summary )


## ------------------------------------------------------------------------
# interpolate both datasets at the time resolution of the neutrophils
dt <- timeStep( TCells )
interpolate.dt <- function( x, dt, how = "spline" ){
  trange <- range( timePoints( wrapTrack( x ) ) )
  tvec <- seq( trange[1], trange[2], by = dt )
  x <- interpolateTrack( x, tvec, how = how )
  return(x)
}

T.corrected <- as.tracks( lapply( TCells, interpolate.dt, dt = dt ) )
T2.corrected <- as.tracks( lapply( T2, interpolate.dt, dt = dt ) )

# Check the effect on the displacement statistics:
T.corr.steps <- subtracks( T.corrected, 1 )
T.corr.disp <- sapply( T.corr.steps, displacement )

T2.corr.steps <- subtracks( T2.corrected, 1 )
T2.corr.disp <- sapply( T2.corr.steps, displacement )

lapply( list( T = T.disp, T2 = T2.disp, 
              T.corr = T.corr.disp, T2.corr = T2.corr.disp ), 
        summary )

