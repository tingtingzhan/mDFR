



#' @title Distribution Free Test Statistic
#' 
#' @description
#' Distribution free test statistic to compare a treatment against a control.
#' 
#' @param x1,x0 \link[base]{numeric} \link[base]{matrix}-es, 
#' treatment \eqn{x_1} and control responses \eqn{x_0}, respectively
#' 
#' @param null.value \link[base]{numeric} scalar \eqn{c}, 
#' as in \eqn{H_0: \mu_1 - c\mu_0 = 0}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @references
#' The distribution free test statistic \eqn{T} is not exactly equation (1) and (2) of Santos' 2015 cell paper \doi{10.3390/cells4010001}.
#' 
#' @keywords internal
#' @importFrom matrixStats rowVars
#' @importFrom stats var
#' @export
santosT <- function(
    x1, x0, 
    null.value = 1,
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
  
  ### Equation (3) of \doi{10.3390/cells4010001} is wrong!!!! 
  ### 'the maximum standard deviation of an n = 6 binary set' # ???
  # sd_pooled <- pmax(sqrt(30)/10, sqrt( ((n1-1L)*vr1 + (n0-1L)*vr0) / (n1+n0-2L)))
  ### reason 1: need to multiply by `null.value`
  ### reason 2: should not assume equal.variance between `x1` and `x0`
  
  delta <- (m1 - null.value * m0)
  
  gw <- Gosset_Welch(v1 = vr1, v0 = vr0, c0 = null.value, n1 = n1, n0 = n0)
  gw_stderr <- attr(gw, which = 'stderr', exact = TRUE)
  
  out <- delta / gw_stderr
  attr(out, which = 'delta') <- delta
  attr(out, which = 'stderr') <- gw_stderr
  attr(out, which = 'df') <- c(gw) # to drop attributes
  return(out)
  
}




#' @title Distribution Free Test Statistic, Difference of Difference
#' 
#' @param u,v two objects, both returned from function [santosT]
#' 
#' @references
#' The distribution free test statistic \eqn{T} is not exactly the same as 
#' Equation (6) and (7) of Santos' 2015 cell paper \doi{10.3390/cells4010001}.
#' 
#' @returns
#' Function [santosT2()] returns a \link[base]{numeric} \link[base]{vector} \eqn{T}.
#' 
#' @export
santosT2 <- function(
    u, v
) {
  d_u <- attr(u, which = 'delta', exact = TRUE)
  d_v <- attr(v, which = 'delta', exact = TRUE)
  sd_u <- attr(u, which = 'stderr', exact = TRUE)
  sd_v <- attr(v, which = 'stderr', exact = TRUE)
  df_u <- attr(u, which = 'df', exact = TRUE)
  df_v <- attr(v, which = 'df', exact = TRUE)
  sd_pooled <- sqrt( (df_u*sd_u^2 + df_v*sd_v^2) / (df_u+df_v) )
  (d_u - d_v) / sd_pooled
}

