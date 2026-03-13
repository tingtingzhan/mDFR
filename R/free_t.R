

#' @title Distribution \linkS4class{free_t}-Statistic
#' 
#' @description
#' Distribution \linkS4class{free_t} statistic to compare a treatment against a control.
#' 
#' @slot .Data \link[base]{numeric} \link[base]{vector}, distribution-free \eqn{t}-statistic
#' 
#' @slot null.value \link[base]{numeric} scalar \eqn{c}

#' @slot delta \link[base]{numeric} \link[base]{vector}
#' 
#' @slot stderr \link[base]{numeric} \link[base]{vector}
#' 
#' @slot df \link[base]{numeric} \link[base]{vector}
#' 
#' @slot data any R object, the input data, for downstream analysis
#' 
#' @references
#' \url{https://tingtingzhan-maxt.netlify.app/appendix/statistic.html#sec-free_t_math}
#' \url{https://tingtingzhan-maxt.netlify.app/s4/free_t#sec-free_t}
#' 
#' @keywords internal
#' @name free_t
#' @aliases free_t-class
#' @export
setClass(Class = 'free_t', contains = 'numeric', slots = c(
  null.value = 'numeric',
  delta = 'numeric',
  stderr = 'numeric',
  df = 'numeric',
  data = 'ANY'
))


setMethod(f = show, signature = 'free_t', definition = \(object) {
  
  'Distribution Free t-Statistic' |>
    message()
  
  object@null.value |> 
    col_blue() |> style_bold() |>
    sprintf(fmt = 'H0: \u03bc1 = c \u03bc0; where c = %s') |>
    message()
  
  object@.Data |>
    print()
  
  'Obtained From' |>
    style_bold() |>
    message()
  object@data |>
    show()
  
})





#' @rdname free_t
#' 
#' @param x see **Usage**
#' 
#' @param pid1 (optional) \link[base]{integer} \link[base]{vector},
#' permuted indices for treatment
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @export
free_t <- function(x, ...) UseMethod(generic = 'free_t')

#' @rdname free_t
#' @export
free_t.ELISpot <- function(x, pid1, ...) {
  
  if (missing(pid1)) {
    return(.free_t(x1 = x@x1, x0 = x@x0, data = x, ...))
  }
  
  z <- cbind(x@x1, x@x0)
  return(.free_t(
    x1 = z[, pid1, drop = FALSE], 
    x0 = z[, -pid1, drop = FALSE], 
    ... # internal use, no need to carry @data
  ))
  
}

