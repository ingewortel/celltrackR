## ----pack, warning = FALSE, message = FALSE------------------------------
library( celltrackR )
library( ggplot2 )

## ----Tdata---------------------------------------------------------------
str( TCells, list.len = 3 )

## ------------------------------------------------------------------------
head( TCells[[1]] )

## ------------------------------------------------------------------------
brownian <- brownianTrack( nsteps = 20, dim = 3 )
str(brownian)

## ------------------------------------------------------------------------
plot( wrapTrack(brownian) )

## ---- fig.width = 7------------------------------------------------------
brownian.tracks <- simulateTracks( 10, brownianTrack( nsteps = 20, dim = 3 ) )
par(mfrow=c(1,2))
plot( brownian.tracks, main = "simulated random walk" )
plot( normalizeTracks( TCells ), main = "real T cell data" )

## ---- fig.width = 7------------------------------------------------------
# get displacement vectors
step.displacements <- t( sapply( subtracks(TCells,1), displacementVector ) )

# get mean and sd of displacement in each dimension
step.means <- apply( step.displacements, 2, mean )
step.sd <- apply( step.displacements, 2, sd )

# simulate brownian motion with the same statistics
brownian.tracks.matched <- simulateTracks( 10, brownianTrack( nsteps = 20, dim = 3,
                                                              mean = step.means,
                                                              sd = step.sd ) )

# compare displacement distributions
data.displacement <- sapply( subtracks( TCells,1), displacement )
matched.displacement <- sapply( subtracks( brownian.tracks.matched, 1 ), displacement )

df <- data.frame( disp = data.displacement,
                  data = "TCells" )
df2 <- data.frame( disp = matched.displacement,
                   data = "model" )
df <- rbind( df, df2 )
ggplot( df, aes( x = data, y = disp ) ) +
  geom_boxplot() +
  theme_classic()

# Plot new simulation versus real data
par(mfrow=c(1,2))
plot( brownian.tracks.matched, main = "simulated random walk" )
plot( normalizeTracks( TCells ), main = "real T cell data" )


## ------------------------------------------------------------------------

# simulate brownian motion with bias
brownian.tracks.bias <- simulateTracks( 10, brownianTrack( nsteps = 20, dim = 3,
                                                              mean = c(1,1,1) ) )

plot( brownian.tracks.bias, main = "biased random walk" )


## ------------------------------------------------------------------------
beauchemin.tracks <- simulateTracks( 10, beaucheminTrack(sim.time=20) )
plot( beauchemin.tracks )

## ------------------------------------------------------------------------
bootstrap.tracks <- simulateTracks( 10, bootstrapTrack( nsteps = 20, TCells ) )
plot( bootstrap.tracks )

## ---- fig.width = 7, warning = FALSE, message = FALSE--------------------
# Simulate more tracks to reduce noice
bootstrap.tracks <- simulateTracks( 100, bootstrapTrack( nsteps = 20, TCells ) )

# Compare step speeds in real data to those in bootstrap data
real.speeds <- sapply( subtracks( TCells,1 ), speed )
bootstrap.speeds <- sapply( subtracks( bootstrap.tracks,1), speed )
dspeed <- data.frame( tracks = c( rep( "data", length( real.speeds ) ),
                                  rep( "bootstrap", length( bootstrap.speeds ) ) ),
                      speed = c( real.speeds, bootstrap.speeds ) )

# Same for turning angles
real.angles <- sapply( subtracks( TCells,2 ), overallAngle, degrees = TRUE )
bootstrap.angles <- sapply( subtracks( bootstrap.tracks,2), overallAngle, degrees = TRUE )
dangle <- data.frame( tracks = c( rep( "data", length( real.angles ) ),
                                  rep( "bootstrap", length( bootstrap.angles ) ) ),
                      angle = c( real.angles, bootstrap.angles ) )

# plot
pspeed <- ggplot( dspeed, aes( x = tracks, y = speed ) ) +
  geom_violin( color = NA, fill = "gray" ) +
  geom_boxplot( width = 0.3 ) +
  theme_classic()

pangle <- ggplot( dangle, aes( x = tracks, y = angle ) ) +
  geom_violin( color = NA, fill = "gray" ) +
  geom_boxplot( width = 0.3 ) +
  theme_classic()

gridExtra::grid.arrange( pspeed, pangle, ncol = 2 )


## ---- fig.width=6--------------------------------------------------------
# Simulate more tracks
brownian.tracks <- simulateTracks( 100, brownianTrack( nsteps = 20, dim = 3,
                                                              mean = step.means,
                                                              sd = step.sd ) )
bootstrap.tracks <- simulateTracks( 100, bootstrapTrack( nsteps = 20, TCells ) )

msd.data <- aggregate( TCells, squareDisplacement, FUN = "mean.se" )
msd.data$data <- "data"
msd.brownian <- aggregate( brownian.tracks, squareDisplacement, FUN = "mean.se" )
msd.brownian$data <- "brownian"
msd.bootstrap <- aggregate( bootstrap.tracks, squareDisplacement, FUN = "mean.se" )
msd.bootstrap$data <-"bootstrap"

msd <- rbind( msd.data, msd.brownian, msd.bootstrap )
ggplot( msd, aes( x = i, y = mean, ymin = lower, ymax = upper, color = data, fill = data ) ) +
  geom_ribbon( color= NA, alpha  = 0.2 ) +
  geom_line() +
  labs( x = "t (steps)",
        y = "square displacement" ) +
  scale_x_continuous(limits= c(0,20) ) +
  theme_bw()

## ---- fig.width=6--------------------------------------------------------
# compute autocorrelation
acor.data <- aggregate( TCells, overallDot, FUN = "mean.se" )
acor.data$data <- "data"
acor.brownian <- aggregate( brownian.tracks, overallDot, FUN = "mean.se" )
acor.brownian$data <- "brownian"
acor.bootstrap <- aggregate( bootstrap.tracks, overallDot, FUN = "mean.se" )
acor.bootstrap$data <-"bootstrap"

acor <- rbind( acor.data, acor.brownian, acor.bootstrap )
ggplot( acor, aes( x = i, y = mean, ymin = lower, ymax = upper, color = data, fill = data ) ) +
  geom_ribbon( color= NA, alpha  = 0.2 ) +
  geom_line() +
  labs( x = "dt (steps)",
        y = "autocovariance" ) +
  scale_x_continuous(limits= c(0,20) ) +
  theme_bw()

