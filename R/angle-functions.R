#' Angle Between Two Vectors
#'
#' a and b can be numeric vectors or matrices in which each row represents a numeric vector.
#' In the last case, one angle is returned for each row.
#'
vecAngle <- function( a, b, degrees = TRUE )
{
  # Checks
  if( class(a) != class(b) ){
    stop( "vecAngle: a and b must both be numeric vectors (of the same length) or matrices (of the same size)")
  }
  if( any( pracma::size(a) != pracma::size(b) ) ){
    stop( "vecAngle: cannot compute angle between vectors of unequal dimensions.")
  }
  if( is.matrix(a) && nrow(a) != nrow(b) ){
    stop( "vecAngle: a and b must have an equal number of rows.")
  }
  
  # Convert vector to matrix
  if( !is.matrix(a) ){
    a <- matrix( a, ncol = length(a) )
    b <- matrix( b, ncol = length(b) )
  }
  
  # Normalize a and b by their length
  a <- a/sqrt(.rowSums(a^2, nrow(a), ncol(a)))
  b <- b/sqrt(.rowSums(b^2, nrow(b), ncol(b)))
  
  # Compute dot product
  rs <- .rowSums(a * b, nrow(a), ncol(a))
  rs[rs > 1] <- 1
  rs[rs < -1] <- -1
  
  # angle is the acos of the dot product of the normalized vectors
  ang <- acos(rs)
  
  if( degrees ){
    return( pracma::rad2deg(ang) )
  } else{
    return(ang)
  }
}


