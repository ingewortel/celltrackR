#' Cluster Tracks
#'
#' Perform a quick clustering visualization of a set of tracks according to a given vector
#' of track measures.
#'
#' @param tracks the tracks that are to be clustered.
#' @param measures a function, or a vector of functions (see \link{TrackMeasures}).
#' Each function is expected to
#' return a single number given a single track.
#' @param scale logical indicating whether the measures values shall be scaled
#' using the function \code{\link[base]{scale}} before the clustering.
#' @param method \code{"hclust"} for hierarchical clustering, or
#'  \code{"kmeans"} for k-means clustering.
#' @param labels optional: a vector of labels of the same length as the track object.
#' These are used to color points in the visualization.
#' @param return.clust logical: return the clustering object instead of only the plot?
#' (defaults to \code{FALSE}).
#' @param ... additional parameters to be passed to the corresponding clustering
#' function: \code{\link[stats]{hclust}} or  \code{\link[stats]{kmeans}}.
#'
#' @return By default, only returns a plot. If \code{return.clust=TRUE}, also returns
#' a clustering object as returned by \code{\link[stats]{hclust}} or  \code{\link[stats]{kmeans}}.
#  See the documentation of those functions for details on the
#'   output object.
#'
#' @details The measures are applied to each of the tracks in the given
#' \emph{tracks} object. According to the resulting values, the tracks are
#' clustered using the chosen clustering method.
#' If \code{scale} is \code{TRUE}, the measure values are scaled to mean value
#' \eqn{0} and standard deviation \eqn{1} (per measure) before the clustering.
#'
#' Method hclust plots a dendrogram of the clustering.
#'
#' Method kmeans plots each computed cluster (x-axis) versus each of the track
#' measures in the \code{measures} vector, producing one panel per measure.
#' If labels are given, points are colored according to their "true" label.
#'
#' @seealso \code{\link{getFeatureMatrix}} to obtain a feature matrix that can be
#' used for manual clustering and plotting, and \code{\link{trackFeatureMap}} to
#' visualize high-dimensional track feature data via dimensionality reduction.
#'
#' @examples
#' ## Cluster tracks according to the mean of their Hust exponents along X and Y
#' ## using hierarchical clustering
#'
#' cells <- c(TCells,Neutrophils)
#' real.celltype <- rep(c("T","N"),c(length(TCells),length(Neutrophils)))
#' ## Prefix each track ID with its cell class to evaluate the clustering visually
#' names(cells) <- paste0(real.celltype,seq_along(cells))
#' clust <- clusterTracks( cells, hurstExponent, method = "hclust",
#'  return.clust = TRUE  )
#'
#' ## How many cells are "correctly" clustered?
#' sum( real.celltype == c("T","N")[cutree(clust,2)] )
#'
#' @export
clusterTracks <- function( tracks, measures, scale = TRUE, labels = NULL, method = "hclust", return.clust = FALSE, ... )
{
  if( length( measures) == 0 ){
    stop( "clusterTracks: no measures given! Please specify at least one.")
  }

  # Get the feature matrix
  values <- getFeatureMatrix( tracks, measures )

  # Scale to mean 0, sd 1 if scale = TRUE
  if (scale) {
    values <- scale(values)
  }

  # Some of the following will adjust the 'par' settings for plotting.
  # Ensure that these will be reset to their default values before the function
  # exits for any reason.
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  # Compute clustering based on "method"
  if( method == "hclust" ){
    clust <- stats::hclust( stats::dist(values), ...)
  } else if( method == "kmeans" ){
    clust <- stats::kmeans( values, ... )
  } else {
    stop( "clusterTracks: unknown method! Please choose either hclust or kmeans." )
  }

  # Make the plot

    # For hierarchical clustering, plot the dendrogram colored by label if labels
    # are given.
  if (method == "hclust") {

    dend <- stats::as.dendrogram(clust)

    # color the leaves of the dendrogram by label
    if( !is.null(labels) ){
      if( !requireNamespace("dendextend", quietly=TRUE ) ){
        stop( "clusterTracks: please install the 'dendextend' package to use this functionality!" )
      }
      colors_to_use <- as.numeric( factor( labels ))
      colors_to_use <- colors_to_use[stats::order.dendrogram(dend)]
      dendextend::labels_colors(dend) <- colors_to_use
    }

    graphics::par( cex = 0.7 )
    graphics::plot( dend )
    if( !is.null(labels) ){
      graphics::legend("topright", legend = unique(labels),
             fill = unique( colors_to_use ),
             border = unique( colors_to_use ), bty = "n")
    }

    # For kmeans clustering, plot each feature in the matrix for each of the k clusters.
    # Color points by label if label is given.
  } else if ( method == "kmeans" ){
    if( is.null(labels) ){
      labs <- 1
    } else {
      labs <- as.numeric( factor(labels) )
    }

    # plotting area
    nc <- 2
    nrow <- ceiling( ncol(values)/nc )

    graphics::par( mfrow=c(nrow, nc ), xpd = TRUE )

    # make plots for each feature
    for( i in 1:ncol(values) ){
      graphics::plot(jitter(clust$cluster), values[,i], type = "p",
           pch = 19, col = labs, xlim = c(0.5, max(clust$cluster)+1),
           main = paste("Feature",i), xlab="cluster", xaxt="n")
      graphics::axis(1, at = seq(1, length(unique( clust$cluster) )))
      if( !is.null(labels) ){
        graphics::legend("topright", legend=unique(labels),
               col=unique(labs), pch = 1, cex=0.8)
      }

    }

  }

  if( return.clust ){
    return(clust)
  }
}

#' Dimensionality Reduction on Track Features
#'
#' Perform a quick dimensionality reduction visualization of a set of tracks according to a given vector
#' of track measures.
#'
#' @param tracks the tracks that are to be clustered.
#' @param measures a function, or a vector of functions (see \link{TrackMeasures}).
#' Each function is expected to
#' return a single number given a single track.
#' @param scale logical indicating whether the measures values shall be scaled
#' using the function \code{\link[base]{scale}} before the mapping is performed.
#' @param method \code{"PCA"} for a
#' scatterplot along principal components, \code{"MDS"} for multidimensional scaling,
#' \code{"UMAP"} for a UMAP. Note that for
#' \code{"UMAP"}, the \code{uwot} package must be installed.
#' @param labels optional: a vector of labels of the same length as the track object.
#' These are used to color points in the visualization.
#' @param return.mapping logical: return the mapping object instead of only the plot?
#' (defaults to \code{FALSE}).
#' @param ... additional parameters to be passed to the corresponding
#' function: \code{\link[stats]{prcomp}} (for \code{method="PCA"}),
#'   \code{\link[stats]{cmdscale}} (for \code{method="MDS"}),
#'   or  \code{\link[uwot]{umap}} (for \code{method="UMAP"}).
#'
#' @return By default, only returns a plot. If \code{return.clust=TRUE}, also returns
#' a clustering object as returned by \code{\link[stats]{hclust}},  \code{\link[stats]{kmeans}},
#'  \code{\link[stats]{prcomp}} (returns \code{$x}), \code{\link[stats]{cmdscale}},
#'   or  \code{\link[uwot]{umap}} (returns \code{$layout}). See the documentation of those functions for details on the
#'   output object.
#'
#' @details The measures are applied to each of the tracks in the given
#' \emph{tracks} object. According to the resulting values, the tracks are
#' mapped to fewer dimensions using the chosen method.
#' If \code{scale} is \code{TRUE}, the measure values are scaled to mean value
#' \eqn{0} and standard deviation \eqn{1} (per measure) before the mapping.
#'
#' The dimensionality reduction methods PCA, MDS, and UMAP each produce a
#' scatterplot of all tracks as points, plotted along the principal component
#' axes generated by the corresponding method.
#'
#' @seealso \code{\link{getFeatureMatrix}} to obtain a feature matrix that can be
#' used for manual clustering and plotting, and \code{\link{clusterTracks}} to
#' perform hierarchical or k-means clustering on a tracks dataset.
#'
#' @examples
#' ## Map tracks according to speed, mean turning angle, straightness, and asphericity
#' ## using multidimensional scaling, and store output.
#'
#' cells <- c(TCells,Neutrophils)
#' real.celltype <- rep(c("T","N"),c(length(TCells),length(Neutrophils)))
#' ## Prefix each track ID with its cell class to evaluate the clustering visually
#' names(cells) <- paste0(real.celltype,seq_along(cells))
#' map <- trackFeatureMap( cells, c(speed,meanTurningAngle,straightness, asphericity),
#'  method = "MDS",  return.mapping = TRUE  )
#'
#' @export
trackFeatureMap <- function( tracks, measures, scale = TRUE, labels = NULL, method = "PCA", return.mapping = FALSE, ... )
{
  if( length( measures) == 0 ){
    stop( "trackFeatureMap: no measures given! Please specify at least one.")
  }

  # Get the feature matrix
  values <- getFeatureMatrix( tracks, measures )

  # Scale to mean 0, sd 1 if scale = TRUE
  if (scale) {
    values <- scale(values)
  }

  # Compute mapping based on "method"
  if( method == "MDS" ){
    clust <- stats::cmdscale( stats::dist(values), ... )
  } else if( method == "PCA" ){
    clust <- stats::prcomp( values, ... )$x
  } else if( method == "UMAP" ){
    if( !requireNamespace("uwot", quietly=TRUE ) ){
      stop( "clusterTracks: please install the 'uwot' package to use this functionality!" )
    }
    clust <- uwot::umap( values, ... )
  } else {
    stop( "trackFeatureMap: unknown method! Please choose from: MDS, PCA, or UMAP." )
  }

  # Make the plot
  # For dimensionality reduction methods, plot points in two dimensions,
  # colored by label if labels are given.
  if( is.element( method, c("MDS","PCA","UMAP") ) ){
    #df <- as.data.frame(clust)
    #colnames(df) <- c("V1","V2")
    if( !is.null(labels) ){
      lab <- as.numeric( factor( labels ) )
    } else {
      lab <- rep(1,length( tracks ) )
    }

    # Before adjusting the 'par' settings for plotting, ensure that they
    # are reset whenever the function exits for any reason.
    oldpar <- graphics::par(no.readonly = TRUE)
	  on.exit(graphics::par(oldpar))

	  # Plot the result
    graphics::par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE )
    graphics::plot( clust, col = lab )
    if( !is.null(labels)){
      graphics::legend("topright", legend=unique(labels), inset=c(-0.1,0),
                       col=unique(lab), pch = 1, cex=0.8)
    }
  }

  if( return.mapping ){
    return(clust)
  }
}


#' Obtaining A Feature Matrix
#'
#' Applies a given vector of track measures directly on a set of tracks, returning
#' output in a matrix with a column for each measure and a row for each track. Can
#' also return a distance matrix, which some clustering methods require.
#'
#' @param tracks the tracks that are to be analyzed.
#' @param measures a function, or a vector of functions (see \link{TrackMeasures}).
#' Each function is expected to return a single number given a single track.
#' @param dist should a distance matrix rather than a feature matrix be returned?
#' @param ... further arguments passed on to "dist"
#'
#' @return A matrix with a row for each track and a column for each measure.
#'
#' @seealso \code{\link{clusterTracks}} for a quick method to compute the feature
#' matrix and a clustering, and \code{\link{trackFeatureMap}} to perform
#' dimensionality reduction methods on a set of track features.
#'
#' @examples
#' ## Get speed, meanTurningAngle, and straightness for T cell tracks
#' fm <- getFeatureMatrix( TCells, c(speed,meanTurningAngle,straightness))
#' str(fm)
#'
#' @export
getFeatureMatrix <- function( tracks, measures, dist = FALSE, ... )
{

  values <- matrix(nrow = length(tracks))
  if (is.function(measures)) {
    measures <- c(measures)
  }
  values <- do.call(cbind, lapply(measures,
                                  function(m) sapply(tracks, m)))

  if( dist ){
    return( dist( values, ... ) )
  } else {
    return(values)
  }

}

