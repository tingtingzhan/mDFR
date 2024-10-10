

#' @title Trim `'elispot'` Object
#' 
#' @param x an `elispot` object
#' 
#' @details
#' Both treatment `$x1` and control `$x0` must contain more than 1 measurement
#' 
#' @returns
#' Function [trim_elispot] returns an `elispot` object.
#' 
#' @export
trim_elispot <- function(x) {
  id <- (rowSums(!is.na(x$x1)) > 1L) & (rowSums(!is.na(x$x0)) > 1L)
  return(x[id, ])
}
