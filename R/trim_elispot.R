

#' @title Trim `'elispot'` Object
#' 
#' @param x an `elispot` object
#' 
#' @details
#' Both treatment `@x1` and control `@x0` must contain more than 1 measurement
#' 
#' @returns
#' Function [trim_elispot()] returns an `elispot` object.
#' 
#' @keywords internal
#' @export
trim_elispot <- function(x) {
  
  id <- (rowSums(!is.na(x@x1)) > 1L) & (rowSums(!is.na(x@x0)) > 1L)
  
  if (all(id)) return(x)
  
  x@design <- x@design[id, ]
  x@x1 <- x@x1[id, , drop = FALSE]
  x@x0 <- x@x0[id, , drop = FALSE]
  return(x)
  
}
