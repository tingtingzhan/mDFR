

#' @title maxT_moodie
#' 
#' @description Resampling-based \eqn{p}-value adjustment
#' 
#' @param data an `elispot` object
#' 
#' @param bootstrap \link[base]{logical} scalar, 
#' whether to use bootstrap samples.
#' Default `FALSE` indicating the use of permuted samples.
#' 
#' @param null.value \link[base]{numeric} scalar \eqn{\mu_0}, as in \eqn{H_0: \bar{X}_1 - \bar{X}_0 = \mu_0}
#' 
# @param log_base \link[base]{numeric} scalar, base of log-transformation.
# Choice of `log_base` does not affect the p-values, only the statistics.
# When `log_base = NULL`, we are testing
# \eqn{H_0: \mu_{e} - \mu_{c} = 0}, i.e., `null.value` is forced to be 0
#' 
#' @param R \link[base]{integer} scalar \eqn{R}, number of bootstrap copies if `bootstrap = TRUE`
#' 
#' @param seed_ (optional) \link[base]{integer} scalar, random seed for `'boot'` resampling method.
#' In Moodie's function `elsdfr2x()`, `seed_ = 9456845L` is used.
#' Howwever, Moodie uses a same seed for all analysis, which is a *bad* practice!
#' The purpose of this parameter is solely to re-create Moodie's results.
#' The end user is advised to leave this parameter missing.
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
#' 
#' @example inst/moodie/maxT_moodie.R
#' 
#' @importFrom methods new
#' @export
maxT_moodie <- function(
    data, bootstrap = FALSE,
    null.value,
    R = 1e3L, seed_, # only if (bootstrap)
    ...
) {
  
  m1 <- rowMeans(data$x1, na.rm = TRUE)
  m0 <- rowMeans(data$x0, na.rm = TRUE)
  t_ <- (m1 - m0) - null.value
  
  if (!bootstrap) { # moodie's [elsdfreq]
    dd <- cbind(data$x1, data$x0)
    ids <- perm_elispot(data)
    prm1. <- lapply(ids, FUN = \(i) rowMeans(dd[, i, drop = FALSE], na.rm = TRUE)) # permutation of treatment
    prm0. <- lapply(ids, FUN = \(i) rowMeans(dd[, -i, drop = FALSE], na.rm = TRUE)) # permutation of control
    prm. <- unlist(prm1., use.names = FALSE) - unlist(prm0., use.names = FALSE)
    dim(prm.) <- c(.row_names_info(data, type = 2L), length(ids))
    #T_ <- prm. # moodie's
    T_ <- prm. - null.value # Tingting's
    # moodie's [elsdfreq] example has `null.value = 0`
    # Tingting thinks moodie's code is wrong
    
  } else { # moodie's [elsdfr2x]
    bmeans <- function(x, R) {
      # moodie's functions [bf] and [bcf], re-written
      # `x` is numeric vector
      x_a <- x[!is.na(x)]
      if (!(nx <- length(x_a))) stop('all-missing not allowed')
      x_b <- sample(x_a, size = nx * R, replace = TRUE)
      rowMeans(array(x_b, dim = c(R, nx)))
    }
    if (!missing(seed_)) set.seed(seed = seed_) 
    bt1 <- t.default(apply(data$x1, MARGIN = 1L, FUN = bmeans, R = R)) # moodie's `bemeans`
    bt0 <- t.default(apply(data$x0, MARGIN = 1L, FUN = bmeans, R = R)) # moodie's `bcmeans`
    T_ <- (bt1 - m1) - (bt0 - m0) # moodie's
    # Tingting does not understand this `T_`
    
  }

  data$x1 <- data$x0 <- NULL
  class(data) <- 'data.frame' # 'elispot' is not specified in slot `@design` of \linkS4class{maxT}
  
  ag0 <- list(...)[c('two.sided')]
  return(do.call(new, args = c(list(
    Class = 'maxT', 
    t. = t_, T. = T_,
    design = data#,
    #name = paste(unlist(tmp[lengths(tmp) == 1L], use.names = FALSE), collapse = '; ')
  ), ag0[lengths(ag0, use.names = FALSE) > 0L])))
  
}






