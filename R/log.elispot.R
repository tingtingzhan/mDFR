

#' @title Several S3 method dispatches for `ELISpot`
#' 
#' @description
#' `ELISpot` dispatches for S3 generics \link[base]{log}.
#' 
#' @param x an `ELISpot` object
#' 
#' @param base \link[base]{numeric} scalar, see \link[base]{log}
#' 
#' @note
#' The `S3` generic function \link[base]{log1p} does not have `base` parameter.
#' 
#' @returns
#' 
#' Function [log.ELISpot()] returns an `ELISpot`.
#' 
#' @keywords internal
#' @name log_ELISpot
#' @export log.ELISpot
#' @export
log.ELISpot <- function(x, base = exp(1)) {
  if (any(x@x1 == 0, x@x0 == 0, na.rm = TRUE)) {
    x@x1 <- x@x1 + 1
    x@x0 <- x@x0 + 1
  }
  x@x0 <- log(x@x0, base = base)
  x@x1 <- log(x@x1, base = base)
  return(x)
}


#' @rdname log_ELISpot
#' @export log1p.ELISpot
#' @export
log1p.ELISpot <- function(x) {
  x@x0 <- log1p(x@x0)
  x@x1 <- log1p(x@x1)
  return(x)
}

#' @rdname log_ELISpot
#' @export log10.ELISpot
#' @export
log10.ELISpot <- function(x) {
  x@x0 <- log10(x@x0)
  x@x1 <- log10(x@x1)
  return(x)
}

#' @rdname log_ELISpot
#' @export log2.ELISpot
#' @export
log2.ELISpot <- function(x) {
  x@x0 <- log2(x@x0)
  x@x1 <- log2(x@x1)
  return(x)
}



