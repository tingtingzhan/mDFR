

#' @title Distribution Free Resampling
#' 
#' @description ..
#' 
#' @param data an `elispot`
#' 
#' @param ... additional parameters for function [split.elispot] and [maxT_moodie]
#' 
#' @references 
#' Original algorithm by Zoe Moodie on March 16, 2010, and
#' her 2006 papers \doi{10.1016/j.jim.2006.07.015},
#' and an empirical study \doi{10.1007/s00262-010-0875-4}.
#' 
#' \url{https://rundfr.fredhutch.org}
#' 
#' @examples 
#' (m0 = moodie2ELISpot(`^a[1-9]$` ~ antigen | day + id, data = moodie, control = 'negctl')) # no log
#' (r1a = DFR(m0, null.value = 0))
#' (r2a = DFR(log(m0, base = 10), null.value = log10(2), resample = 'boot'))
#' 
#' @export
DFR <- function(data, ...) {
  ds <- split.elispot(data, ...)
  ret0 <- lapply(ds, FUN = function(i, ...) {
    maxT_moodie(i, ...)
  }, ...)
  ret <- do.call(rbind.data.frame, args = c(ret0, list(make.row.names = FALSE)))
  # stopifnot(identical(class(ret), c('DFR', 'data.frame')))
  return(ret)
}

#' @export
print.DFR <- function(x, digits = 3L, row.names = FALSE, ...) {
  x0 <- x
  x0$y1 <- x0$y0 <- NULL
  print.data.frame(x0, digits = digits, row.names = FALSE, ...)
}


#' @title maxT_moodie
#' 
#' @description Resampling-based \eqn{p}-value adjustment
#' 
#' @param data an `elispot`
#' 
#' @param resample \link[base]{character} scalar, `'combn'` (default) or `'boot'`
#' 
#' @param null.value \link[base]{numeric} scalar \eqn{\mu_0}, as in \eqn{H_0: \bar{x}_1 - \bar{x}_0 = \mu_0}
#' 
# @param log_base \link[base]{numeric} scalar, base of log-transformation.
# Choice of `log_base` does not affect the p-values, only the statistics.
# When `log_base = NULL`, we are testing
# \eqn{H_0: \mu_{e} - \mu_{c} = 0}, i.e., `null.value` is forced to be 0
#' 
#' @param R \link[base]{integer} scalar \eqn{R}, number of bootstrap copies if `resample = 'boot'`
#' 
#' @param seed_ (optional) \link[base]{integer} scalar, random seed for `'boot'` resampling method.
#' In Moodie's function `elsdfr2x()`, `seed_ = 9456845L` is used.
#' Howwever, Moodie uses a same seed for all analysis, which is a *bad* practice!
#' The purpose of this parameter is solely to re-create Moodie's results.
#' The end user is advised to leave this parameter missing.
#' 
#' @param two.sided see [maxT_]
#' 
#' @param ... potential parameters
#' 
#' @details ..
#' 
#' @references 
#' ISBN: 978-0-471-55761-6
#' 
#' \url{https://rundfr.fredhutch.org}
#' 
#' @export
maxT_moodie <- function(
    data, resample = c('combn', 'boot'),
    null.value,
    R = 1e3L, seed_, # only for `resample == 'boot'`
    two.sided = FALSE,
    ...
) {
  
  m1 <- rowMeans(data$y1, na.rm = TRUE)
  m0 <- rowMeans(data$y0, na.rm = TRUE)
  t_ <- (m1 - m0) - null.value
  
  switch(match.arg(resample), combn = { # moodie's [elsdfreq]
    dd <- cbind(data$y1, data$y0)
    ids <- combn_elispot(data)
    # permutations
    prm1 <- do.call(cbind, args = lapply(ids, FUN = function(i) rowMeans(dd[, i, drop = FALSE], na.rm = TRUE)))
    prm0 <- do.call(cbind, args = lapply(ids, FUN = function(i) rowMeans(dd[, -i, drop = FALSE], na.rm = TRUE)))
    resamp <- prm1 - prm0
    
  }, boot = { # moodie's [elsdfr2x]
    bmeans <- function(x, R) {
      # moodie's functions [bf] and [bcf], re-written in a beautiful way
      # `x` is numeric vector
      x0 <- x[!is.na(x)]
      if (!(nx <- length(x0))) stop('all-missing not allowed')
      x1 <- sample(x0, size = nx * R, replace = TRUE)
      rowMeans(array(x1, dim = c(R, nx)))
    }
    if (!missing(seed_)) set.seed(seed = seed_) 
    bt1 <- t.default(apply(data$y1, MARGIN = 1L, FUN = bmeans, R = R)) # moodie's `bemeans`
    bt0 <- t.default(apply(data$y0, MARGIN = 1L, FUN = bmeans, R = R)) # moodie's `bcmeans`
    resamp <- (bt1 - m1) - (bt0 - m0)
  })
  
  ret <- data.frame(
    data,
    tstat = t_,
    adjp = maxT_(t. = t_, T. = resamp, two.sided = two.sided)
  )
  class(ret) <- c('DFR', class(ret))
  return(ret)

}

# ?utils::combn is not S3 generic
#' @importFrom utils combn
combn_elispot <- function(data) {
  # no randomness
  n1 <- dim(data$y1)[2L]
  n0 <- dim(data$y0)[2L]
  # if ((n1 > 5L) && (n0 > 5L)) message('combn slow for large number of replicates!')
  ret <- combn(n1 + n0, m = n1, simplify = FALSE)
  if (length(ret) < 20L) stop('Too few combn replicates to use this method (< 20).')
  return(ret)
}

