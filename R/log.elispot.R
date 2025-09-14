

#' @title Several S3 method dispatches for `elispot`
#' 
#' @description
#' `elispot` dispatches for S3 generics \link[base]{log}.
#' 
#' @param x an `elispot` object
#' 
#' @param base \link[base]{numeric} scalar, see \link[base]{log}
#' 
#' @note
#' The `S3` generic function \link[base]{log1p} does not have `base` parameter.
#' 
#' @returns
#' 
#' Function [log.elispot()] returns an `elispot`.
#' 
#' @keywords internal
#' @name log_elispot
#' @export log.elispot
#' @export
log.elispot <- function(x, base = exp(1)) {
  if (any(x@x1 == 0, x@x0 == 0, na.rm = TRUE)) {
    x@x1 <- x@x1 + 1
    x@x0 <- x@x0 + 1
  }
  x@x0 <- log(x@x0, base = base)
  x@x1 <- log(x@x1, base = base)
  return(x)
}


#' @rdname log_elispot
#' @export log1p.elispot
#' @export
log1p.elispot <- function(x) {
  x@x0 <- log1p(x@x0)
  x@x1 <- log1p(x@x1)
  return(x)
}

#' @rdname log_elispot
#' @export log10.elispot
#' @export
log10.elispot <- function(x) {
  x@x0 <- log10(x@x0)
  x@x1 <- log10(x@x1)
  return(x)
}

#' @rdname log_elispot
#' @export log2.elispot
#' @export
log2.elispot <- function(x) {
  x@x0 <- log2(x@x0)
  x@x1 <- log2(x@x1)
  return(x)
}



