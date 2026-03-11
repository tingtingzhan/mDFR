

#' @title Create \linkS4class{free_t} Object
#' 
#' @description
#' Low-level utility function to create \linkS4class{free_t} and \linkS4class{free_d} objects.
#' 
#' @param x1 \link[base]{numeric} \link[base]{matrix},
#' treatment responses \eqn{x_1}
#' 
#' @param x0 \link[base]{numeric} \link[base]{matrix},
#' control responses \eqn{x_0}
#' 
#' @param null.value \link[base]{numeric} scalar, the constant \eqn{c},
#' default value is 1 
#' 
#' @param s4 \link[base]{logical} scalar, whether to return an `S4` object,
#' default value is `TRUE`. Set as `FALSE` for internal use to save compute time.
#' 
#' @param data (optional) an R object, 
#' the first argument of function [free_d] or [free_t],
#' for downstream analysis
#' 
#' @param ... additional parameters, currently of no use
#' 
#' @references
#' \url{https://tingtingzhan-maxt.netlify.app/S4/free_t.html#sec-free_t_create}
#' 
#' @name free_internal
#' @importFrom matrixStats rowVars
#' @importFrom stats var
.free_t <- function(
    x1, x0, null.value = 1,
    data,
    s4 = TRUE,
    ...
) {
  
  d1 <- dim(x1)
  m1 <- x1 |>
    .rowMeans(m = d1[1L], n = d1[2L], na.rm = TRUE)
  n1 <- (!is.na(x1)) |>
    .rowSums(m = d1[1L], n = d1[2L], na.rm = TRUE)
  vr1 <- x1 |>
    rowVars(na.rm = TRUE)
  # ?matrixStats::rowVars is the fastest solution I know of!!!
  
  d0 <- dim(x0)
  m0 <- x0 |>
    .rowMeans(m = d0[1L], n = d0[2L], na.rm = TRUE)
  n0 <- (!is.na(x0)) |>
    .rowSums(m = d0[1L], n = d0[2L], na.rm = TRUE)
  vr0 <- x0 |>
    rowVars(na.rm = TRUE)
  
  delta <- (m1 - null.value * m0)
  
  gw <- Gosset_Welch(v1 = vr1, v0 = vr0, c0 = null.value, n1 = n1, n0 = n0)
  gw_stderr <- attr(gw, which = 'stderr', exact = TRUE)
  
  out <- delta / gw_stderr
  attr(out, which = 'null.value') <- null.value
  attr(out, which = 'delta') <- delta
  attr(out, which = 'stderr') <- gw_stderr
  attr(out, which = 'df') <- c(gw) # to drop attributes
  
  attr(out, which = 'data') <- if (!missing(data)) data # else NULL
  
  if (s4) return(new(Class = 'free_t', out)) # `slot`s are `attr`s !
  
  return(out)
  
}



#' @rdname free_internal
.free_d <- function(
    x1, x0, null.value = 0,
    data,
    s4 = TRUE,
    ...
) {
  
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


