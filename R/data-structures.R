#' Tracks Objects
#'
#' The function \code{tracks} is used to create tracks objects. \code{as.tracks} coerces
#' its argument to a tracks object, and \code{is.tracks} tests for tracks objects.
#' \code{c} can be used to combine (concatenate) tracks objects.
#'
#' @param x an object to be coerced or tested.
#'
#' @param ... for \code{tracks}, numeric matrices or objects that can be coerced to
#'  numeric matrices. Each
#'  matrix contains the data of one track. The first column is the time, and the remaining
#'  columns define a spatial position. Every given matrix has to contain the same number
#'  of columns, and at least two columns are necessary.
#'
#'  For \code{c}, tracks objects to be combined.
#'
#'  For \code{as.tracks}, further arguments passed to methods (currently not used).
#'
#' @details Tracks objects are lists of matrices. Each matrix contains at least two
#' columns; the first column is time, and the remaining columns are a spatial coordinate.
#' The following naming conventions are used (and enforced by \code{tracks}): The time
#' column has the name `t`, and spatial coordinate columns have names `x`,`y`,`z` if there
#' are three or less coordinates, and `x1`,...,`xk` if there are \eqn{k \ge 4}
#' coordinates. All tracks in an object must have the same number of dimensions. The
#' positions in a track are expected to be sorted by time (and the constructor
#' \code{tracks} enforces this).
#'
#' @examples
#' ## A single 1D track
#' x <- tracks( matrix(c(0, 8,
#' 10, 9,
#' 20, 7,
#' 30, 7,
#' 40, 6,
#' 50, 5), ncol=2, byrow=TRUE ) )
#'
#' ## Three 3D tracks
#' x2 <- tracks( rbind(
#'  c(0,5,0), c(1,5,3), c(2,1,3), c(3,5,6) ),
#'  rbind( c(0,1,1),c(1,1,4),c(2,5,4),c(3,5,1),c(4,-3,1) ),
#'  rbind( c(0,7,0),c(1,7,2),c(2,7,4),c(3,7,7) ) )
#'
#' @name tracks
NULL


#' Convert Tracks to Data Frame
#'
#' Converts tracks from the list-of-matrices format, which is good
#' for efficient processing and therefore the default in this package, to a
#' single dataframe which is convenient for plotting or saving the data.
#'
#' @param x the \code{tracks} object to be coerced to a data frame.
#'
#' @param row.names NULL or a character vector giving row names for the
#'  data frame.  Missing values are not allowed.
#' @param optional logical. Required for S3 consistency, but
#' has no effect: column names are always assigned to the resulting
#'  data frame regardless of the setting of this option.
#' @param include.timepoint.column logical. If set to \code{TRUE}, then the resulting
#'  dataframe will contain a column that consecutively numbers the positions according
#'  to their time. Note that this information is anyway implicitly present in the time
#'  information.
#' @param idsAsFactors logical. If \code{TRUE}, then the id column of the resulting
#'  dataframe will be a factor column, otherwise a characeter column.
#' @param ... further arguments to be passed from or to other methods.
#'
#' @return A single data frame containing all individual tracks from the input with a
#' prepended column named "id" containing each track's identifier in `x`.
#'
#' @examples
#' ## Display overall average position of the T cell data
#' colMeans( as.data.frame( TCells )[-c(1,2)] )
#' @export
as.data.frame.tracks <- function(x, row.names = NULL, optional = FALSE,
                                 include.timepoint.column=FALSE, idsAsFactors = TRUE, ...) {
  ids <- rep(names(x), sapply(x,nrow))
  if( include.timepoint.column ){
    timepoint <- stats::ave( ids, ids, FUN=seq_along )
    r <- data.frame(id=ids, timepoint=timepoint, do.call(
      rbind.data.frame, x ), stringsAsFactors=idsAsFactors )
  } else {
    r <- data.frame(id=ids, do.call(
      rbind.data.frame, x ), stringsAsFactors=idsAsFactors )
  }
  if( !is.null(row.names) && length(row.names)==nrow(r) ){
    rownames(r) <- row.names
  } else {
    rownames(r) <- NULL
  }
  return(r)
}

#' Convert from Data Frame to Tracks
#'
#' Get cell tracks from a data.frame. Data are expected to be organized as
#' follows.
#' One column contains a track identifier, which can be numeric or a string, and
#' determines which points belong to the same track.
#' Another column is expected to contain a time index or a time period (e.g. number of
#' seconds elapsed since the beginning of the track, or since the beginning of the
#' experiment). Input of dates is not (yet) supported, as absolute time information is
#' frequently not available.
#' Further columns contain the spatial coordinates. If there are three or less spatial
#' coordinates, their names will by "x", "y", and "z"
#' (depending on whether the tracks are 1D, 2D or 3D). If there are four or more spatial
#' coordinates, their names will be "x1", "x2", and so on.
#' The names or indices of these columns in the data.frame are given using the
#' corresponding parameters (see below). Names and indices can be mixed, e.g. you can
#' specify \code{id.column="Parent"} and \code{pos.columns=1:3}
#'
#' @inheritParams read.tracks.csv
#' @param x the data frame to be coerced to a \code{tracks} object.
#'
#' @return A \code{tracks} object.
#'
#' @export
as.tracks.data.frame <- function(x, id.column=1, time.column=2,
                                 pos.columns=c(3,4,5), scale.t=1,
                                 scale.pos=1, ...) {
  if( ncol(x) < length(pos.columns) + 1) {
    stop("Data frame does not contain enough columns! (Perhaps you need to specify 'sep')")
  }
  if( length(pos.columns) < 1 ){
    stop("At least one position column needs to be specified!")
  }

  # Special case: if columns are in form e.g. c("x",NA) then we
  # read all columns from the beginning to the end.
  if( length(pos.columns) == 2 && is.na(pos.columns[2]) ){
    cx <- match( pos.columns[1], colnames(x) )
    if( is.na(cx) && is.numeric(pos.columns[1]) ){
      cx <- pos.columns[1]
    }
    pos.columns <- seq( cx, ncol(x) )
  }
  cx <- as.character(c(id.column,time.column,pos.columns))
  cxc <- match( cx, colnames(x) )
  cxi <- match( cx, seq_len(ncol(x)) )

  cxi[is.na(cxi)] <- cxc[is.na(cxi)]

  if( any(is.na(cxi)) ){
    stop("Column(s) not found: ",
         paste(cx[is.na(cxi)],collapse=","))
  }
  r <- x[,as.integer(cxi)]
  if( ncol(r) <= 5 ){
    colnames(r) <- c("id","t",c("x","y","z")[seq_along(pos.columns)])
  } else {
    colnames(r) <- c("id","t",paste0("x",seq_along(pos.columns)))
  }
  if( scale.t != 1 ){
    r[,"t"] <- scale.t*r[,"t"]
  }
  if( any( scale.pos != 1 ) ){
    r[,-c(1,2)] <- scale.pos*r[,-c(1,2)]
  }
  sort.tracks(as.tracks.list(split.data.frame(as.matrix(r[,-1]), r[,1])))
}

#' @export
`[.tracks` <- function(x,y) as.tracks(as.list(x)[y])

#' @rdname tracks
#' @return A \code{tracks} object.
#' @export
as.tracks <- function(x, ...)
  UseMethod("as.tracks")

#' @rdname tracks
#' @export
as.tracks.list <- function(x, ...)
  structure(x, class="tracks")


#' Convert from Tracks to List
#'
#' Coerces a \code{tracks} object to a list.
#'
#' @param x the \code{tracks} object to be coerced to a list.
#' @param ... further arguments to be passed from or to other methods.
#' @return A generic list of single tracks, where each track is a matrix with 
#' \code{t/delta.t} rows and 4 columns. This looks a lot like a tracks object,
#' except that its class is not "tracks" anymore.	
#' @export
as.list.tracks <- function(x, ...)
  structure( x, class="list" )


#' @rdname tracks
#' @export
is.tracks <- function(x)
  inherits(x, "tracks")


#' @rdname tracks
#' @export
c.tracks <- function(...) {
  args <- lapply(list(...), as.list)
  as.tracks(do.call(c, args))
}

#' Create Track Object from Single Track
#'
#' Makes a \code{tracks} object containing the given track.
#'
#' @param x the input track.
#'
#' @return A list of class \code{tracks} containing only the input track \code{x}, which
#' is assigned the name "1".
#' @export
wrapTrack <- function(x) {
  return(as.tracks(list("1" = x)))
}

#' @rdname tracks
#' @export
tracks <- function( ... ){
  tlist <- lapply( list(...), data.matrix )
  if( !all( sapply( tlist , is.numeric ) ) ){
    stop("Track time/position information must be numeric!")
  }
  tdims <- sapply( tlist, ncol )
  if( any( tdims < 2 ) ){
    stop("At least time and one spatial dimension are needed!")
  }
  if( length(tlist) == 0 ){
    return( as.tracks(tlist) )
  }
  if( length(tlist) > 0 &&
      diff( range( tdims ) ) > 0 ){
    stop("All tracks must have the same number of dimensions!")
  }
  r <- tlist[[1]]
  if( ncol(r) <= 5 ){
    cls <- c("t",c("x","y","z")[seq_len(ncol(r)-1)])
    for( i in seq_along(tlist) ){
      colnames(tlist[[i]]) <- cls
    }
  } else {
    cls <- c("t",paste0("x",seq_len(ncol(r)-1)))
    for( i in seq_along(tlist) ){
      colnames(tlist[[i]]) <- cls
    }
  }
  if( is.null(names(tlist)) ){
    names(tlist) <- seq_along(tlist)
  }
  sort(as.tracks(tlist))
}
