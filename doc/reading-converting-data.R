## ----pack, warning = FALSE, message = FALSE------------------------------
library( celltrackR )
library( ggplot2 )

## ------------------------------------------------------------------------
d <- read.table( system.file("extdata", "t-cells.txt", package="celltrackR" ) )
str(d)
head(d)

## ------------------------------------------------------------------------
t <- read.tracks.csv( system.file("extdata", "t-cells.txt", package="celltrackR" ), 
		      header = FALSE, 
                      id.column = 2, time.column = 3, pos.columns = 4:6 )
plot(t)

## ------------------------------------------------------------------------
# Structure of the TCells object
str( t, list.len = 3 )

# This object is both a list and a "tracks" object
is.list( t )
is.tracks( t )

# The first element is the track of the first cell in the data:
head( t[[1]] )

## ------------------------------------------------------------------------
# Get the first track
t1 <- t[[1]]
str(t1)

# This is no longer a tracks object, but a matrix
is.tracks( t1 )
is.matrix( t1 )

## ----fig.width =7--------------------------------------------------------
par( mfrow=c(1,2))
plot( t1, main = "Plotting matrix directly" )
plot( wrapTrack( t1 ), main = "After using wrapTrack()" )

## ------------------------------------------------------------------------
# Get the first track
t1b <- t[1]
str(t1b)

# This remains a track object
is.tracks( t1b )

## ------------------------------------------------------------------------
# Get the first and the third track
t13 <- t[c(1,3)]
str(t13)

## ------------------------------------------------------------------------
# Get tracks with ids 1 and 3
t13b <- t[c("1","3")]
str(t13b)

## ------------------------------------------------------------------------
speeds <- sapply( t, speed )
speeds

## ------------------------------------------------------------------------
# Function to remove all data after given timepoint
# x must be a single track matrix, which is what this function will
# receive from lapply
removeAfterT <- function( x, time.cutoff ){
  
  # Filter out later timepoints
  x2 <- x[ x[,"t"] <= time.cutoff, ]
  
  # Return the new matrix, or NULL if there are no timepoints before the cutoff
  if( nrow(x2) == 0 ){
    return(NULL)
  } else {
    return(x2)
  }
}

# Call function on each track using lapply
filtered.t <- lapply( t, function(x) removeAfterT( x, 200 ) )

# Remove any tracks where NULL was returned
filtered.t <- filtered.t[ !sapply( filtered.t, is.null )]

## ------------------------------------------------------------------------
str(filtered.t, list.len = 3 )
is.list( filtered.t )
is.tracks( filtered.t )

## ---- fig.width = 7------------------------------------------------------
filtered.t <- as.tracks( filtered.t )
is.tracks( filtered.t )

par(mfrow=c(1,2))
plot( t, main = "Unfiltered data")
plot( filtered.t, main = "Filtered on timepoints < 200" )

## ----fig.width = 7-------------------------------------------------------
# The filtering function must return TRUE or FALSE for each track given to it
my.filter <- function(x){
  return( nrow(x) > 15 )
}

# Filter with this function using filterTracks
long.tracks <- filterTracks( my.filter, t )

# Plot the result
par(mfrow=c(1,2))
plot( t, main = "All tracks")
plot( long.tracks, main = "Long tracks only" )

## ----fig.width = 7-------------------------------------------------------

# Filter with this function using filterTracks
median.speed <- median( sapply( t, speed ) )
fast.tracks <- selectTracks( t, speed, median.speed, Inf )

# Plot the result
par(mfrow=c(1,2))
plot( t, main = "All tracks")
plot( fast.tracks, main = "Fastest half" )

## ---- fig.width = 7------------------------------------------------------
# Lower resolution
lower.res <- subsample( t, k = 2 )

# Plot the result
par(mfrow=c(1,2))
plot( t, main = "Original data")
plot( lower.res, main = "Lower resolution" )

## ------------------------------------------------------------------------
subtrack.nsteps <- 2
t.2steps <- subtracks( t, subtrack.nsteps )
str( t.2steps, list.len = 3 )

## ------------------------------------------------------------------------
# Last step of the first subtrack and first step of the second are equal
t.2steps[c(1,2)]

## ------------------------------------------------------------------------
t.2steps.b <- subtracks( t, subtrack.nsteps, overlap = 0 )

# No longer any overlap
t.2steps.b[c(1,2)]

## ------------------------------------------------------------------------
t.prefixes <- prefixes( t, subtrack.nsteps )

# these subtracks come from different cells
t.prefixes[c(1,2)]


## ------------------------------------------------------------------------
# Check which timepoints occur in the dataset
tp <- timePoints(t)
tp

# Extract all subtracks starting from the third timepoint
t.sbytime <- subtracksByTime( t, tp[3], subtrack.nsteps )

t.sbytime[c(1,2)]

## ------------------------------------------------------------------------
# Original tracks object
str( t, list.len = 3 )

# Converted to dataframe
t.df <- as.data.frame(t)
str( t.df )

# Converted to list (note class at the bottom)
t.list <- as.list(t)
str( t.list, list.len = 3 )

# Convert list back to tracks
str( as.tracks( t.list ), list.len = 3 )

# Convert dataframe to tracks
str( as.tracks( t.df ), list.len = 3 )


