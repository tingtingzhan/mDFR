

#' @title Distribution Free Test Statistic \linkS4class{free_t}
#' 
#' @description
#' Distribution free test statistic to compare a treatment against a control.
#' 
#' @slot .Data ..
#' 
#' @slot null.value \link[base]{numeric} scalar \eqn{c}, 
#' as in \eqn{H_0: \mu_1 - c\mu_0 = 0}

#' @slot delta ..
#' 
#' @slot stderr ..
#' 
#' @slot df ..
#' 
#' @slot data ..
#' 
#' @references
#' The distribution free test statistic \eqn{T} is not exactly equation (1) and (2) of Santos' 2015 cell paper \doi{10.3390/cells4010001}.
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


#' @title Distribution Free Test Statistic, Difference of Difference \linkS4class{free_t_diff}
#' 
#' @slot .Data \link[base]{numeric} \link[base]{vector}
#' 
#' @slot e1,e2 \linkS4class{free_t}
#' 
#' @name free_t_diff
#' @aliases free_t_diff-class
#' @export
setClass(Class = 'free_t_diff', contains = 'numeric', slots = c(
  e1 = 'free_t',
  e2 = 'free_t'
))



#' @rdname free_t
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
free_t <- function(x, ...) UseMethod(generic = 'free_t')

#' @rdname free_t
#' @export free_t.ELISpot
#' @export
free_t.ELISpot <- function(x, id1, ...) {
  
  if (missing(id1)) {
    return(free_t.matrix(x = x@x1, x0 = x@x0, data = x, ...))
  }
  
  z <- cbind(x@x1, x@x0)
  return(free_t.matrix(
    x = z[, id1, drop = FALSE], 
    x0 = z[, -id1, drop = FALSE], 
    ...
  ))
  
}


#' @rdname free_t
#' @importFrom matrixStats rowVars
#' @importFrom stats var
#' @export free_t.matrix
#' @export
free_t.matrix <- function(
    x, x0, 
    data,
    null.value = 1,
    s4 = TRUE,
    ...
) {
  
  x1 <- x; x <- NULL
  
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
  
  ### Equation (3) of \doi{10.3390/cells4010001} is wrong!!!! 
  ### 'the maximum standard deviation of an n = 6 binary set' # ???
  # sd_pooled <- pmax(sqrt(30)/10, sqrt( ((n1-1L)*vr1 + (n0-1L)*vr0) / (n1+n0-2L)))
  ### reason 1: need to multiply by `null.value`
  ### reason 2: should not assume equal.variance between `x1` and `x0`
  
  delta <- (m1 - null.value * m0)
  
  gw <- Gosset_Welch(v1 = vr1, v0 = vr0, c0 = null.value, n1 = n1, n0 = n0)
  gw_stderr <- attr(gw, which = 'stderr', exact = TRUE)
  
  out <- delta / gw_stderr
  attr(out, which = 'null.value') <- null.value
  attr(out, which = 'delta') <- delta
  attr(out, which = 'stderr') <- gw_stderr
  attr(out, which = 'df') <- c(gw) # to drop attributes
  attr(out, which = 'x1') <- x1
  attr(out, which = 'x0') <- x0
  attr(out, which = 'data') <- if (!missing(data)) data # else NULL
  if (s4) return(new(Class = 'free_t', out)) # `slot`s are `attr`s !
  return(out)
  
}



#' @title Distribution Free Test Statistic, Difference of Difference
#' 
#' @param e1,e2 \linkS4class{free_t}
#' 
#' @references
#' The distribution free test statistic \eqn{T} is not exactly the same as 
#' Equation (6) and (7) of Santos' 2015 cell paper \doi{10.3390/cells4010001}.
#' 
#' @keywords internal
#' @export
setMethod(f = '-', signature = c(e1 = 'free_t', e2 = 'free_t'), definition = \(e1, e2) {
  
  if (!all.equal.numeric(e1@null.value, e2@null.value)) stop('`@null.value` must be the same')
  
  # `@` much faster than ?base::attr !!!
  d1 <- e1@delta
  d2 <- e2@delta
  sd1 <- e1@stderr
  sd2 <- e2@stderr
  df1 <- e1@df
  df2 <- e2@df
  sd_pooled <- sqrt( (df1*sd1^2 + df2*sd2^2) / (df1+df2) )
  
  new(
    Class = 'free_t_diff',
    (d1 - d2) / sd_pooled,
    e1 = e1, e2 = e2
  )

})

