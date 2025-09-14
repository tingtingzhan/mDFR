
#' @title [split.elispot()]
#' 
#' @param x an \linkS4class{elispot}
#' 
#' @param f ..
#' 
#' @param drop \link[base]{logical} scalar, must be `FALSE`
#' 
#' @keywords internal
#' @export split.elispot
#' @export
split.elispot <- function(x, f, drop = TRUE, ...) {
  
  if (!drop) stop('`drop` must be TRUE')
  
  # most simple solution..
  
  z <- x@design
  z$x1 <- x@x1
  z$x0 <- x@x0

  foo <- \(y) {
    y0 <- y
    y0$x1 <- y0$x0 <- NULL
    new(Class = 'elispot', design = y0, x1 = y$x1, x0 = y$x0)
  }
  
  
  z |>
    split.data.frame(f = f, drop = drop) |>
    lapply(FUN = foo)

}