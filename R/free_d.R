

#' @title Distribution \linkS4class{free_d}-Statistic
#' 
#' @description
#' Distribution \linkS4class{free_d} statistic to compare a treatment against a control.
#' 
#' @slot .Data \link[base]{numeric} \link[base]{vector}, distribution-free \eqn{d}-statistic
#' 
#' @slot null.value \link[base]{numeric} scalar \eqn{c}
#' 
#' @slot data any R object, the input data, for downstream analysis
#' 
#' @references
#' \url{https://tingtingzhan-maxt.netlify.app/appendix/statistic.html#sec-free_d_math}
#' \url{https://tingtingzhan-maxt.netlify.app/S4/free.html#sec-free_d}
#' 
#' @keywords internal
#' @name free_d
#' @aliases free_d-class
#' @export
setClass(Class = 'free_d', contains = 'numeric', slots = c(
  null.value = 'numeric',
  data = 'ANY'
))




#' @rdname free_d
#' 
#' @param x see **Usage**
#' 
#' @param pid1 (optional) \link[base]{integer} \link[base]{vector},
#' permuted indices for treatment (via the function [permID()])
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @export
free_d <- function(x, ...) UseMethod(generic = 'free_d')

#' @rdname free_d
#' @export
free_d.ELISpot <- function(x, pid1, ...) {
  
  if (missing(pid1)) {
    return(.free_d(x1 = x@x1, x0 = x@x0, data = x, ...))
  }
  
  z <- cbind(x@x1, x@x0)
  return(.free_d(
    x1 = z[, pid1, drop = FALSE], 
    x0 = z[, -pid1, drop = FALSE], 
    ... # internal use, no need to carry @data
  ))
  
}

