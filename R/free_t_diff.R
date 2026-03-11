
#' @title Distribution Free \eqn{t}-Statistic, Difference of Difference 
#' 
#' @description
#' The `S4` class \linkS4class{free_t_diff}.
#' 
#' @slot .Data \link[base]{numeric} \link[base]{vector}
#' 
#' @slot e1,e2 \linkS4class{free_t} objects
#' 
#' @references 
#' \url{https://tingtingzhan-maxt.netlify.app/S4/free_t.html#sec-free_t_diff}
#' 
#' @keywords internal
#' @export
setClass(Class = 'free_t_diff', contains = 'numeric', slots = c(
  e1 = 'free_t',
  e2 = 'free_t'
))




#' @title Difference of \linkS4class{free_t} Objects
#' 
#' @description
#' Difference of two \linkS4class{free_t} objects that have a common slot `@null.value`.
#' 
#' @param e1,e2 two \linkS4class{free_t} objects that have a common slot `@null.value`.
#' 
#' @references 
#' \url{https://tingtingzhan-maxt.netlify.app/S4/free.html#sec-free_t_diff}
#' 
#' @keywords internal
#' @export
setMethod(f = '-', signature = c(e1 = 'free_t', e2 = 'free_t'), definition = \(e1, e2) {
  
  if (!all.equal.numeric(e1@null.value, e2@null.value)) stop('`@null.value` must be the same')
  
  # `@` much faster than ?base::attr !!!
  sd1 <- e1@stderr
  sd2 <- e2@stderr
  df1 <- e1@df
  df2 <- e2@df
  sd_pooled <- sqrt( (df1*sd1^2 + df2*sd2^2) / (df1+df2) )
  
  new(
    Class = 'free_t_diff',
    (e1@delta - e2@delta) / sd_pooled,
    e1 = e1, e2 = e2
  )
  
})
