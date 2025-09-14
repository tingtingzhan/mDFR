



#' @title Distribution Free Test Statistic
#' 
#' @description
#' Distribution free test statistic to compare a treatment against a control.
#' 
#' @param x1,x0 \link[base]{numeric} \link[base]{matrix}-es, 
#' treatment \eqn{x_1} and control responses \eqn{x_0}, respectively
#' 
#' @param null.value \link[base]{numeric} scalar \eqn{c}, 
#' as in \eqn{H_0: \bar{X}_1 - c\bar{X}_0 = 0}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @returns
#' Function [santosT] returns a \link[base]{numeric} \link[base]{vector} of \eqn{\bar{x}_1 - c\bar{x}_0}, 
#' i.e., the numerator of distribution free test statistic \eqn{t},
#' with \link[base]{attributes}
#' \describe{
#' \item{`attr(.,'df')`}{degree of freedom \eqn{\text{df}}, to be used in function [santosT2];}
#' \item{`attr(.,'stderr')`}{standard error \eqn{s_{\bar\Delta}}, i.e., the denominator of distribution free test statistic \eqn{t};}
#' \item{`attr(.,'stderr2')`}{standard deviation \eqn{s^2_{\bar\Delta}}, included due to the extreme time-sensitivity of function [santosT2].}
#' }
#' 
#' @references
#' The distribution free test statistic \eqn{T} is not exactly equation (1) and (2) of Santos' 2015 cell paper \doi{10.3390/cells4010001}.
#' 
#' @details
#' Distribution free test statistic
#' \deqn{
#' T = \dfrac{\bar{X}_1 - c\bar{X}_0}{s_{\bar\Delta}}
#' }
#' Let \eqn{n_1}, \eqn{n_0}, \eqn{s^2_1} and \eqn{s^2_0} be sample sizes and sample standard deviations of treatment and control, respectively.
#' The standard deviation 
#' \eqn{s^2_{\bar\Delta}= s^2_1/n_1 + c^2s^2_0/n_0},
#' has degree of freedom based on Welch–Satterthwaite equation (see more from help files of function \link[DanielBiostatistics10th]{Gosset_Welch}),
#' \deqn{\text{df} = \dfrac{\left(\dfrac{s^2_1}{n_1}+\dfrac{c^2s^2_0}{n_0}\right)^2}{\dfrac{(s^2_1/n_1)^2}{n_1-1}+\dfrac{(c^2s^2_0/n_0)^2}{n_0-1}}}
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
  
  tmp <- Gosset_Welch(v1 = vr1, v0 = vr0, c0 = null.value, n1 = n1, n0 = n0)
  ret <- (m1 - null.value * m0)
  attr(ret, which = 'df') <- c(tmp) # to drop attributes
  attr(ret, which = 'stderr') <- attr(tmp, which = 'stderr', exact = TRUE)
  attr(ret, which = 'stderr2') <- attr(tmp, which = 'stderr2', exact = TRUE)
  return(ret)
  
}




#' @title Distribution Free Test Statistic, Difference of Difference
#' 
#' @param u,v two objects, both returned from function [santosT]
#' 
#' @details
#' At time point \eqn{u} (see more at [santosT]), distribution free test statistic
#' \deqn{
#' T_u = \dfrac{\bar{X}_{1,u} - c_u\bar{X}_{0,u}}{s_{\bar\Delta,u}}
#' }
#' The standard deviation 
#' \eqn{s^2_{\bar\Delta,u}=s^2_{1,u}/n_{1,u} + c_u^2s^2_{0,u}/n_{0,u}}
#' has degree of freedom based on Welch–Satterthwaite equation
#' \deqn{\text{df}_u = \dfrac{\left(\dfrac{s^2_{1,u}}{n_{1,u}}+\dfrac{c_u^2s^2_{0,u}}{n_{0,u}}\right)^2}{\dfrac{(s^2_{1,u}/n_{1,u})^2}{n_{1,u}-1}+\dfrac{(c_u^2s^2_{0,u}/n_{0,u})^2}{n_{0,u}-1}}}
#' Similarly, at time point \eqn{v}, we have distribution free test statistic \eqn{T_v},
#' standard deviation \eqn{s^2_{\bar\Delta,v}} and degree of freedom \eqn{\text{df}_v}.
#' 
#' To compare time point \eqn{u} vs. \eqn{v}, we have distribution free test statistic
#' \deqn{
#' T = \dfrac{\left(\bar{X}_{1,u} - c_u\bar{X}_{0,u}\right) - \left(\bar{X}_{1,v} - c_v\bar{X}_{0,v}\right)}{s_{\bar\Delta}}
#' }
#' where 
#' \deqn{
#' s^2_{\bar\Delta} = \dfrac{\text{df}_u\cdot s^2_{\bar\Delta,u} + \text{df}_v\cdot s^2_{\bar\Delta,v}}{\text{df}_u + \text{df}_v}
#' }
#' 
#' @note
#' In practice we often let \eqn{c=c_u=c_v} for more intuitive interpretation.
#' 
#' @references
#' The distribution free test statistic \eqn{T} is not exactly the same as 
#' Equation (6) and (7) of Santos' 2015 cell paper \doi{10.3390/cells4010001}.
#' 
#' @returns
#' Function [santosT2] returns a \link[base]{numeric} \link[base]{vector} \eqn{T}.
#' 
#' @export
santosT2 <- function(
    u, v
) {
  sd2_u <- attr(u, which = 'stderr2', exact = TRUE)
  sd2_v <- attr(v, which = 'stderr2', exact = TRUE)
  df_u <- attr(u, which = 'df', exact = TRUE)
  df_v <- attr(v, which = 'df', exact = TRUE)
  sd_pooled <- sqrt( (df_u*sd2_u + df_v*sd2_v) / (df_u+df_v) )
  (u - v) / sd_pooled
}

