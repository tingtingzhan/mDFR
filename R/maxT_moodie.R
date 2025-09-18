

#' @title maxT_moodie
#' 
#' @description Resampling-based \eqn{p}-value adjustment
#' 
#' @param data an \linkS4class{ELISpot} object
#' 
#' @param null.value \link[base]{numeric} scalar \eqn{\mu_0}, as in \eqn{H_0: \log\mu_1 - \log\mu_0 = \mu_0}
#' 
#' @param R \link[base]{integer} scalar \eqn{R}, number of bootstrap copies if `bootstrap = TRUE`
#' 
#' @param two.sided \link[base]{logical} scalar, see \linkS4class{maxT}
#' 
#' @param ... potential parameters
#' 
#' @details ..
#' 
#' @references 
#' 
#' Original algorithm by Zoe Moodie on March 16, 2010 \url{https://rundfr.fredhutch.org}, 
#' and
#' her 2006 papers \doi{10.1016/j.jim.2006.07.015},
#' and an empirical study \doi{10.1007/s00262-010-0875-4}.
#' 
#' @keywords internal
#' @export
maxT_moodie <- function(
    data,
    null.value,
    R = 1e3L,
    two.sided = TRUE, 
    ...
) {
  
  m1 <- rowMeans(data@x1, na.rm = TRUE)
  m0 <- rowMeans(data@x0, na.rm = TRUE)
  t_ <- (m1 - m0) - null.value
  
  { # moodie's [elsdfr2x]; bootstrap
    bmeans <- \(x, R) {
      # moodie's functions [bf] and [bcf], re-written
      # `x` is numeric vector
      x_a <- x[!is.na(x)]
      if (!(nx <- length(x_a))) stop('all-missing not allowed')
      sample(x_a, size = nx * R, replace = TRUE) |> 
        array(dim = c(R, nx)) |>
        rowMeans()
    }
    bt1 <- data@x1 |>
      apply(MARGIN = 1L, FUN = bmeans, R = R) |>
      t.default() # moodie's `bemeans`
    bt0 <- data@x0 |>
      apply(MARGIN = 1L, FUN = bmeans, R = R) |>
      t.default() # moodie's `bcmeans`
    T_ <- (bt1 - m1) - (bt0 - m0) # moodie's
    # Tingting does not understand this `T_`
    
  }

  new(
    Class = 'maxT', 
    t. = t_, T. = T_,
    design = data@design,
    two.sided = two.sided
  )
}




# call moodie's statistic `free_d`

