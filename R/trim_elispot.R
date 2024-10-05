

#' @title Trim `'elispot'` Object
#' 
#' @param x an `elispot` object
#' 
#' @details
#' Both treatment `$y1` and control `$y0` must contain more than 1 measurement
#' 
#' @returns
#' Function [trim_elispot] returns an `elispot` object.
#' 
#' @export
trim_elispot <- function(x) {
  id <- (rowSums(!is.na(x$y1)) > 1L) & (rowSums(!is.na(x$y0)) > 1L)
  return(x[id, ])
}