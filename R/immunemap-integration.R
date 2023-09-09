#' Read tracks from ImmuneMap
#'
#' Reads tracks from \url{https://immunemap.org} for import into celltrackR. 
#' This produces both tracks object(s) and a dataframe with metadata.
#'
#' @param url of the json file to download from immunemap; this should be the url to the
#'		video metadata without the "/tracks" suffix. With this method, the metadata will
#'		be used to automatically scale time to seconds and coordinates to microns if
#'		\code{scale.auto=TRUE}.
#' @param tracks.url optional: alternatively, provide directly the url of the tracks (ending with "/tracks"),
#' 		or an url of a local json file with tracks. With this method, scales must be set
#'		manually. If not specified, it is assumed that adding the suffix "/tracks" to the
#'		supplied \code{url} will provide the track data.
#' @param input the output of \code{parse.immap.json} serves as input for \code{get.immap.tracks}
#' @param keep.id logical: keep track ids from immunemap? If false, new unique ids are
#'  generated. Defaults to \code{TRUE}. If there are no ids in the input json, a warning
#'  will be returned; this can be suppressed by setting keep.id = \code{FALSE}.
#' @param scale.auto logical: if \code{TRUE} (the default), scales will be set automatically using
#'		the metadata found in \code{url}. This works only if the \code{url} is given, not
#'		if only \code{tracks.url} is supplied.
#' @param scale.t optional: multiply timepoints with constant factor to rescale time.
#' 	By default, immunemap returns time in # frames.
#' @param scale.pos optional: multiply coordinates with constant factor to rescale lengths.
#'  By default, immunemap measures coordinates in pixels.
#' @param warn.scaling logical: if \code{scale.t} and \code{scale.pos} are not set,
#' 	warn the user that units are pixels and #frames instead of microns and min/sec. 
#'  Defaults to \code{TRUE}.
#' @param simplify.2D logical: if \code{TRUE} (default), automatically project to 2D when the
#'  z-coordinate has only one value.
#' @param split.celltypes logical: if \code{TRUE} (default = \code{FALSE}), return not one
#'  tracks object but a list of tracks objects for each celltype in the data (as 
#'  determined from the metadata in the immunemap json).
#' @param warn.celltypes logical: if \code{TRUE} (default), warn when the user is either 
#'  trying to return a single tracks object while the metadata indicates there are 
#'  multiple celltypes in the data, or when the user is trying to set \code{split.celltypes = TRUE} 
#'  when there is only one celltype present.
#' @param ... additional parameters to be passed to \code{\link{get.immap.metadata}}.
#'
#' @return \code{read.immap.json}  returns a list with:
#' \item{tracks}{either a single tracks object or a named list of tracks objects per cell type (if \code{split.celltypes = TRUE}}
#' \item{metadata}{a dataframe with metadata for all the track.ids; this is read from the immunemap json file.}
#' 
#' \code{parse.immap.json} simply returns the R list generated from the input json file.
#' 
#' \code{get.immap.tracks} returns a single tracks object.
#'
#' @details \code{read.immap.json} internally uses \code{parse.immap.json} to parse the json file,
#' \code{get.immap.tracks} to extract the tracks, and  \code{\link{get.immap.metadata}}
#' to read the metadata. 
#' 
#' @note This functionality requires the rjson package to be installed.
#'
#' @seealso \code{\link{get.immap.metadata}}.
#'
#' @examples
#' \dontrun{
#' ## Read tracks from immunemap online, using the video info for automatic scaling
#' tr <- read.immap.json( url = "https://api.immunemap.org/video/14" )
#' 
#' ## Read tracks and rescale time (.5min/frame) and coordinates (2microns/pixel)
#' tracksUrl <- "https://api.immunemap.org/video/14/tracks"
#' tr <- read.immap.json( tracks.url = tracksUrl, scale.auto = FALSE, scale.t = .5, scale.pos = 2 )
#' }
#' 
#' ## Read tracks from a file 
#' # tr <- read.immap.json( tracks.url = "my-file.json", warn.scaling = FALSE )
#'
#' @name ReadImmuneMap
#'
#' @export
read.immap.json <- function( url, tracks.url = NULL, keep.id = TRUE, scale.auto = TRUE, scale.t = NULL, scale.pos = NULL, warn.scaling = TRUE, simplify.2D = TRUE, warn.celltypes = TRUE, split.celltypes = FALSE, ... ){

	if( is.null( tracks.url ) ){ tracks.url <- paste0( url, "/tracks" ) }

	# Read json from file or url; error if not json
	input <- parse.immap.json( tracks.url )
	
	# Check format of the input list.
	.check.immap.json( input )

	# Now we can read the tracks:
	if( scale.auto ){
		video.info <- parse.immap.json( url )
		fps <- video.info$size$fps
		voxel.size <- video.info$size$spacing
		tracks <- get.immap.tracks( input, keep.id = keep.id, scale.t = (1/fps), scale.pos = voxel.size, simplify.2D = simplify.2D )
	} else {
		tracks <- get.immap.tracks( input, keep.id = keep.id, scale.t = scale.t, scale.pos = scale.pos, warn.scaling = warn.scaling, simplify.2D = simplify.2D )
	}	
	
	# And the metadata, parsing arguments ... on 
	meta.df <- get.immap.metadata( input, ... )
	if( !is.element( "id", colnames(meta.df) ) ){
		# if no ids, add just a number and put this column first. 
		meta.df$id <- 1:nrow(meta.df)
		meta.df[,c(ncol(meta.df),1:(ncol(meta.df)-1))]
	}
	meta.df$id <- as.character( meta.df$id )
	
	# Depending on settings split tracks by celltypes :
	celltypes <- unique( meta.df$cellTypeName )
	if( split.celltypes ){
		if( length( celltypes ) == 1 ){
			if( warn.celltypes ) warning( "Cannot split.celltypes when there is only one; returning a single tracks object." )
		} else {
			ids.by.type <- lapply( celltypes, function(x) as.character( meta.df$id[ meta.df$cellTypeName == x] ) )
			tracks.by.type <- lapply( ids.by.type, function(x) tracks[x] )
			names( tracks.by.type ) <- celltypes
			tracks <- tracks.by.type
		}
	} else {
		if( length( celltypes ) > 1 ){
			if( warn.celltypes ) warning( "Returning a single tracks object but there are multiple cellTypeNames; please check metadata!" )
		}
	}
	
	# Return a list with 'tracks' and 'metadata'.
	return(list( tracks = tracks, metadata = meta.df ))

}

.check.immap.json <- function( json.input ){
	# Input must be a list, elements must correspond with tracks (check.immap.single)
	if( !is.list(json.input) ){
		stop( "Error reading from immunemap. Expecting a list of tracks, please check the format." )
	}
	elements.check <- sapply( json.input, .check.immap.single, error = FALSE )
	if( !all( elements.check ) ){
		stop( "Error reading json from ImmuneMap: each track in the json file should be an object that contains a key 'points'. Please check json format." )
	}
	# Check also the points
	points.len <- sapply( json.input, function(x) length(x$points) > 1 )
	if( !any( points.len  ) ) stop( "Error in reading track data: your tracks all contain either a single coordinate or are empty." )
	if( !all( points.len ) ) {
		num.missing <- sum( !points.len )
		message( paste0( "...skipping ",num.missing," tracks with no coordinates or only a single coordinate" ))
	}
	
	points.check <- sapply( json.input, function(x) .check.immap.points( x$points, error = FALSE ) )
	if( !all( points.check ) ){
		stop( "Error reading json from ImmuneMap: the 'points' key in the json object should contain an array of all numeric arrays of length 4. Some elements do not fulfill this criterion; please check format." )
	}
}

.check.immap.single <- function( track.json, error = TRUE ){
	
	# minimum requirement is that the list contains 'points'.
	if( !is.element( "points", names( track.json ) ) ){
		if(error) stop( "Error in reading json from ImmuneMap: each track in the json should be an object that must at least contains a key 'points'. Please check json format." )
		return(FALSE)
	}
	return(TRUE)
}

.check.immap.points <- function( pts, error = TRUE ){

	# 'points' should be a list and should contain at least one element.
	stop.msg <- "Error in reading json from ImmuneMap: 'points' should contain an array of all numeric arrays of length 4. Your 'points' don't fit this format - please check."
	if( !is.list( pts ) ){
		if( error ) stop( stop.msg )
		return (FALSE)
	} 
	#if( length( pts ) == 0 ){
	#	if( error ) stop( stop.msg )
	#	return( FALSE )
	#}

	# 'points' should be a list of numeric vectors of each length 4.
	# (if there is no z coordinate, it still contains a 1 on that position.)
	stop.msg <-  "Error in reading json from ImmuneMap: the 'points' key in the json object should contain an array of all numeric arrays of length 4. Some elements do not fulfill this criterion; please check format."
		
	check.numerics <- sapply( pts, is.numeric )
	if( all( check.numerics ) ){
		check.length <- sapply( pts, length )
		if( any( check.length != 4 ) ){
			if (error) stop( stop.msg )
			return(FALSE)
		}
	}
	
	if( any( !check.numerics ) ){
		if (error ) stop( stop.msg )
		return(FALSE)
	}
	return(TRUE)
	
}

#' @rdname ReadImmuneMap
#' @export
parse.immap.json <- function( url ){
	if( !requireNamespace("rjson", quietly=TRUE ) ){
      stop( "Trying to read tracks from ImmuneMap: please install the 'rjson' package to use this functionality!" )
    }

	input <- tryCatch( rjson::fromJSON( file = url ), 
		error = function(cond){ 
			message(paste("Error reading url/file:", url))
			message("Are you sure this is a json file? Here's the original error message:")
           stop( cond )
	} )
	return(input)
}

#' @rdname ReadImmuneMap
#' @export
get.immap.tracks <- function( input, keep.id = TRUE, scale.t = NULL, scale.pos = NULL, warn.scaling = TRUE, simplify.2D = TRUE){

	
	tracks <- lapply( input, .read.immap.single, keep.id = keep.id, scale.t = scale.t, scale.pos = scale.pos, warn.scaling = warn.scaling )
	tracks <- as.tracks( unlist( tracks, recursive = FALSE ) )
	
	# If simplify.2D, remove last coordinate if it is the same everywhere.
	if( simplify.2D ){
		z.coord <- unlist( sapply( tracks, function(x) x[,4] ) )
		if( length( unique( z.coord == 1 ) ) ) tracks <- projectDimensions( tracks )
	}
	
	return(tracks)
}

.read.immap.single <- function( track.json, keep.id = TRUE, scale.t = NULL, scale.pos = NULL, warn.scaling = TRUE ){
	
	# check that this is correct format for a track, return error otherwise.
	.check.immap.single( track.json ) 
	
	# check format of the 'points', return error if problem
	.check.immap.points( track.json$points )
	
	# if points are empty or contain a single coordinate; return no track (NULL)
	if( length( track.json$points) < 2 ){
		return( NULL )
	}
	
	# Read points
	tx <- matrix( unlist( track.json$points ), ncol = 4, byrow = TRUE )
	colnames(tx) <- c("t","x","y","z")
	
	# Default units are 'pixels' and 'steps'. To get to microns and a time unit, 
	# scale.t and scale.pos must be supplied. If not supplied, give a warning 
	# unless this is turned off.
	if( warn.scaling ){
		if( is.null( scale.pos ) ){
			warning( "In reading tracks from ImmuneMap: spatial scale of data unnkown, using pixels. Set parameter 'scale.pos' to supply the spatial resolution, or turn off this warning using 'warn.scaling=FALSE'." )
		}
		if( is.null( scale.t ) ){
			warning( "In reading tracks from ImmuneMap: temporal scale of data unnkown, using frames. Set parameter 'scale.t' to supply the time step between frames, or turn off this warning using 'warn.scaling=FALSE'." )
		}
	}
	# if they are supplied, apply the scaling.
	if( !is.null( scale.t ) ) tx[,1] <- tx[,1]*scale.t
	if( !is.null( scale.pos ) ) tx[,2:4] <- t( t(tx[,2:4])*scale.pos )
	
	# if keep.id = TRUE the track should contain an id.
	id <- NULL
	if( keep.id ){
		if( !is.element( "id", names(track.json) ) ){
			warning( "In reading tracks from ImmuneMap json: keep.id is set to TRUE but the track contains no id. Returning a track without id. To avoid this message, set keep.id = FALSE." )
		} else {
			id <- as.character( track.json$id )
		}
	}	
	
	# Make track
	out <- wrapTrack( tx )
	if( !is.null(id) ) names(out) <- id
	return(out)
	
}


#' Get  Track Metadata from ImmuneMap 
#'
#' Get metadata from tracks obtained from \url{https://immunemap.org} and import into celltrackR. 
#'
#' @param input a parsed json file obtained with \code{\link{parse.immap.json}}
#' @param warn.exclude logical: if \code{TRUE} (default), warn when key-value pairs in the json 
#'  (other than those in exclude.names) are being ignored while parsing immunemap json.
#' @param exclude.names if the json contains keys with these names, they are ignored when reading
#'  the metadata. 
#'
#' @return a dataframe with metadata. This function currently only handles metadata with a single
#' 	value for each track and ignores others (with a warning when \code{warn.exclude=TRUE}).
#' 	column names in the dataframe correspond to the keys in the original json, and values to
#'  the values for each track. 
#'
#' @examples
#' ## Read tracks from immunemap online
#' input <- parse.immap.json( url = "https://api.immunemap.org/video/14/tracks" )
#' meta.df <- get.immap.metadata( input )
#' 
#' ## Repeat but ignore also the 'color' column:
#' exclude <-  c("points", "cellTypeObject","date", "color")
#' meta.df <- get.immap.metadata( input, exclude.names = exclude )
#'
#' @export
get.immap.metadata <- function( input, warn.exclude = TRUE, exclude.names = c("points", "cellTypeObject", "date" ) ){
	
	# Also read the metadata.
	meta.df <- lapply( input, .read.immap.single.metadata, exclude.names = exclude.names, warn.exclude = warn.exclude )
	meta.df <- do.call( rbind, meta.df )
	return(meta.df)
}

.read.immap.single.metadata <- function( track.json, warn.exclude = TRUE, exclude.names = c("points", "cellTypeObject", "date" ) ){
	
	# Take all names except exclude.names from list
	nm <- names( track.json )[ !is.element( names(track.json), exclude.names )  ]
	metadata <- track.json[ nm ]
	
	# Check if these are all single values, then convert to dataframe.
	if( all( sapply( metadata, length ) == 1 ) ){
		metadata <- as.data.frame( metadata )
	} 
	# Otherwise remove anything that doesn't fit and give a warning.
	else {
		ignored <- names( metadata )[ sapply( metadata, length ) != 1  ]
		if( length(ignored) != 0 ){
			metadata <- metadata[ !is.element( names( metadata ), ignored ) ]
			if ( warn.exclude ) warning( paste0("Ignoring tag(s): [ ", ignored,  " ] in track metadata." ) )
		}
	}
	
	# Add the date back if it's in there
	if( nrow( metadata ) == 0 ){
		metadata <- data.frame( date = NA )
	} else{
		metadata$date <- NA
	}
	
	if( is.element( "date", names( track.json ) ) ){
		metadata$date <- paste( track.json$date$date, track.json$date$timezone )
	}
	
	return(metadata)
	
}

    

