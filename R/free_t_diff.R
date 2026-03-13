
#' @title Distribution Free \eqn{t}-Statistic, Difference of Difference 
#' 
#' @description
#' The `S4` class \linkS4class{free_t_diff}.
#' 
#' @slot .Data \link[base]{numeric} \link[base]{vector}
#' 
#' @slot e1,e2 \linkS4class{free_t} objects
#' 
#' @slot design \link[base]{data.frame}, study design showing the difference-of-difference
#' 
#' @references 
#' \url{https://tingtingzhan-maxt.netlify.app/S4/free_t.html#sec-free_t_diff}
#' 
#' @keywords internal
#' @export
setClass(Class = 'free_t_diff', contains = 'numeric', slots = c(
  e1 = 'free_t',
  e2 = 'free_t',
  design = 'data.frame'
))


setMethod(f = show, signature = 'free_t_diff', definition = \(object) {
  
  'Distribution Free t-Statistic for Difference-of-Difference' |>
    message()
  
  object@e1@null.value |> 
    col_blue() |> style_bold() |>
    sprintf(fmt = 'H0: (\u03bc1(\u03d5) - c \u03bc0(\u03d5)) - (\u03bc1(\u03c8) - c \u03bc0(\u03c8)) = 0; where c = %s') |>
    message()
  
  object@.Data |>
    print()
  
  'Obtained From Two' |>
    style_bold() |>
    message()
  object@e1@data |>
    show()
  object@e2@data |>
    show()
  
})



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
  
  if (nrow(e1@data@design) != nrow(e2@data@design)) {
    # timepoint1 and timepoint2 may not have same subjects!!
    #d_1 <- x@e1@data@design
    #d_0 <- x@e2@data@design
    #tmp <- mapply(FUN = intersect, d_1, d_0)
    #d_share <- as.data.frame.list(tmp[lengths(tmp) > 0L])
    stop('um, tzh needs to think about how to write this beautifully..')
    #data1 <- merge.data.frame(d_share, x@e1@data, by = names(d_share), all.x = TRUE)
    #data0 <- merge.data.frame(d_share, x@e2@data, by = names(d_share), all.x = TRUE)
  }
  
  if (!identical(names(e1@data@design), names(e2@data@design))) stop('`@design` at two time points must be the same')
  
  # `@` much faster than ?base::attr !!!
  sd1 <- e1@stderr
  sd2 <- e2@stderr
  df1 <- e1@df
  df2 <- e2@df
  sd_pooled <- sqrt( (df1*sd1^2 + df2*sd2^2) / (df1+df2) )
  
  mapply(
    FUN = \(c1, c0) { # operation per-*c*olumn
      if (anyNA(c1) || anyNA(c0)) stop('does not allow NA in experiment design')
      if (all(c1 == c0)) return(c1)
      return(paste(c1, c0, sep = ' vs. '))
    }, 
    c1 = e1@data@design, 
    c0 = e2@data@design, 
    SIMPLIFY = FALSE
  ) |>
    as.data.frame.list() |>
    new(
      Class = 'free_t_diff',
      (e1@delta - e2@delta) / sd_pooled,
      e1 = e1, e2 = e2,
      design = _
    )

})


#' @export
labels.free_t_diff <- function(object, ...) {
  object@design |>
    labels_design()
}

