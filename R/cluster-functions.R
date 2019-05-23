#' Cluster Tracks
#' 
#' Perform a hierarchical clustering of a set of tracks according to a given vector
#' of track measures.
#' 
#' @param tracks the tracks that are to be clustered.
#' @param measures a function, or a vector of functions (see \link{TrackMeasures}). 
#' Each function is expected to 
#' return a single number given a single track.
#' @param scale logical indicating whether the measures values shall be scaled
#' using the function \code{\link[base]{scale}} before the clustering. 
#' @param ... additional parameters to be passed to \code{\link[stats]{hclust}}.
#'
#' @return An object of class *hclust*, see \code{\link[stats]{hclust}}.
#' 
#' @details The measures are applied to each of the tracks in the given
#' \emph{tracks} object. According to the resulting values, the tracks are 
#' clustered using a hierarchical clustering (see \code{\link[stats]{hclust}}).
#' If \code{scale} is \code{TRUE}, the measure values are scaled to mean value 
#' \eqn{0} and standard deviation \eqn{1} (per measure) before the clustering.
#'
#' Compute clustering
#' supported methods: hclust, kmeans, cmdscale, pca, umap 
#' @examples 
#' ## Cluster tracks according to the mean of their Hust exponents along X and Y
#'
#' cells <- c(TCells,Neutrophils)
#' real.celltype <- rep(c("T","N"),c(length(TCells),length(Neutrophils)))
#' ## Prefix each track ID with its cell class to evaluate the clustering visually
#' names(cells) <- paste0(real.celltype,seq_along(cells))
#' clust <- clusterTracks( cells, hurstExponent )
#' plot( clust )
#' ## How many cells are "correctly" clustered?
#' sum( real.celltype == c("T","N")[cutree(clust,2)] )
#' 
clusterTracks <- function( tracks, measures, scale = TRUE, labels = NULL, method = "hclust", return.clust = FALSE, ... ) 
{
  # Get the feature matrix
  values <- getFeatureMatrix( tracks, measures )
  
  # Scale to mean 0, sd 1 if scale = TRUE
  if (scale) {
    values <- scale(values)
  }
  
  # Compute clustering based on "method"
  if( method == "hclust" ){
    clust <- hclust( dist(values), ...)
  } else if( method == "MDS" ){
    clust <- cmdscale( dist(values), ... )
  } else if( method == "kmeans" ){
    clust <- kmeans( values, ... )
  } else if( method == "PCA" ){
    clust <- prcomp( values, ... )$x
  } else if( method == "UMAP" ){
    clust <- umap::umap( values )$layout
  } else {
    stop( "clusterTracks: unknown method! Please choose from: hclust, MDS, kmeans, PCA, or UMAP." )
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
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE )
    plot( clust, col = lab )
    if( !is.null(labels)){ 
      legend("topright", legend=unique(labels), inset=c(-0.1,0),
             col=unique(lab), pch = 1, cex=0.8)
    }
    par( mar=c(5.1, 4.1, 4.1, 2.1), xpd = FALSE )
    
    # For hierarchical clustering, plot the dendrogram colored by label if labels
    # are given.
  } else if (method == "hclust") {
    
    dend <- as.dendrogram(clust)
    
    # color the leaves of the dendrogram by label
    if( !is.null(labels) ){
      colors_to_use <- as.numeric( factor( labels ))
      colors_to_use <- colors_to_use[order.dendrogram(dend)]
      dendextend::labels_colors(dend) <- colors_to_use
    }
    
    plot( dend )
    if( !is.null(labels) ){
      legend("topright", legend = unique(labels),
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
    par( mfrow=c(nrow, nc ), xpd = TRUE )
    
    # make plots for each feature
    for( i in 1:ncol(values) ){
      plot(jitter(clust$cluster), values[,i], type = "p",
           pch = 19, col = labs,
           main = paste("Feature",i), xlab="cluster", xaxt="n")
      axis(1, at = seq(1, length(unique( clust$cluster) )))
      if( !is.null(labels) ){
        legend("topright", legend=unique(labels),
               col=unique(labs), inset=c(-0.2,0), pch = 1, cex=0.8)        
      }
      
    }
    
    # reset par
    par( mfrow=c(1,1), xpd = FALSE)
  }
  
  if( return.clust ){
    return(clust)
  }
}

#' Obtaining A Feature Matrix
#' 
#' Applies a given vector of track measures directly on a set of tracks, returning
#' output in a matrix with a column for each measure and a row for each track.
#' 
#' @param tracks the tracks that are to be analyzed.
#' @param measures a function, or a vector of functions (see \link{TrackMeasures}). 
#' Each function is expected to return a single number given a single track.
#' 
#' @return A matrix with a row for each track and a column for each measure.
#' 
#' @examples 
#' ## Get speed, meanTurningAngle, and straightness for T cell tracks
#' fm <- getFeatureMatrix( TCells, c(speed,meanTurningAngle,straightness))
#' str(fm)
#' 
#'
getFeatureMatrix <- function( tracks, measures )
{
  
  values <- matrix(nrow = length(tracks))
  if (is.function(measures)) {
    measures <- c(measures)
  }
  values <- do.call(cbind, lapply(measures, 
                                  function(m) sapply(tracks, m)))
  
  return(values)
}

