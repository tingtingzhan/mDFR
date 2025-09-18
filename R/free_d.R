

#' @title Distribution \linkS4class{free_d}-Statistic
#' 
#' @description
#' Distribution \linkS4class{free_d} statistic to compare a treatment against a control.
#' 
#' @slot .Data ..
#' 
#' @slot null.value \link[base]{numeric} scalar \eqn{c}, 
#' as in \eqn{H_0: \mu_1 - c\mu_0 = 0}
#' 
# @slot delta ..
# 
# @slot stderr ..
# 
# @slot df ..
#' 
#' @slot data ..
#' 
#' @keywords internal
#' @name free_d
#' @aliases free_d-class
#' @export
setClass(Class = 'free_d', contains = 'numeric', slots = c(
  null.value = 'numeric',
  #delta = 'numeric',
  #stderr = 'numeric',
  #df = 'numeric',
  data = 'ANY'
))




#' @rdname free_d
#' 
#' @param x see **Usage**
#' 
#' @param x0 \link[base]{numeric} \link[base]{matrix},
#' control responses \eqn{x_0}
#' 
#' @param id1 (optional) \link[base]{integer} \link[base]{vector},
#' permuted indices for treatment
#' 
#' @param null.value see **Slots**
#' 
#' @param s4 \link[base]{logical} scalar
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @export
free_d <- function(x, ...) UseMethod(generic = 'free_d')

#' @rdname free_d
#' @export free_d.ELISpot
#' @export
free_d.ELISpot <- function(x, id1, ...) {
  
  if (missing(id1)) {
    return(free_d.matrix(x = x@x1, x0 = x@x0, data = x, ...))
  }
  
  z <- cbind(x@x1, x@x0)
  return(free_d.matrix(
    x = z[, id1, drop = FALSE], 
    x0 = z[, -id1, drop = FALSE], 
    ... # internal use, no need to carry @data
  ))
  
}


#' @rdname free_d
#' @importFrom matrixStats rowVars
#' @importFrom stats var
#' @export free_d.matrix
#' @export
free_d.matrix <- function(
    x, x0, 
    data,
    null.value = 0,
    s4 = TRUE,
    ...
) {
  
  x1 <- x; x <- NULL
  
  d1 <- dim(x1)
  m1 <- x1 |>
    .rowMeans(m = d1[1L], n = d1[2L], na.rm = TRUE)
  
  d0 <- dim(x0)
  m0 <- x0 |>
    .rowMeans(m = d0[1L], n = d0[2L], na.rm = TRUE)
  
  out <- (m1 - m0) - null.value
  attr(out, which = 'null.value') <- null.value
  attr(out, which = 'data') <- if (!missing(data)) data # else NULL
  if (s4) return(new(Class = 'free_d', out)) # `slot`s are `attr`s !
  return(out)
  
}



