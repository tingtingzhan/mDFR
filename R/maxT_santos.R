
#' @title Modified Distribution Free Resampling
#' 
#' @description ..
#' 
#' @param data an `elispot`
#' 
#' @param ... additional parameters, such as `null.value` in function [santosT],
#' and `two.sided` for function [maxT]
#' 
#' @examples
#' # paper page 5, 'The ELISPOT spot results were analyzed for statistical significance at 
#' #  .. a given time point by comparing the number of spots for each of the tested peptides 
#' #  .. (triplicates) with the background controls (no peptide added, six replicates) 
#' #  .. at that time point.'
#' # paper page 8, layout of Figure 2B 
#' # these are reasons why we use `f = ~ Hr + antigen`
#' st1 = split(santos1, f = ~ Hr + antigen)
#' st2 = split(santos2, f = ~ Hr + antigen)
#' maxT_santos(st1[[1L]], null.value = 1)
#' maxT_santos(st1[[1L]], null.value = 2)
#' lapply(st1, FUN = maxT_santos, null.value = 1)
#' lapply(st1, FUN = maxT_santos, null.value = 2)
#' lapply(st2, FUN = maxT_santos, null.value = 1)
#' lapply(st2, FUN = maxT_santos, null.value = 2)
#' 
#' @references 
#' 
#' Radleigh Santos's 2015 paper \doi{10.3390/cells4010001}
#' 
#' @export
maxT_santos <- function(data, ...) { 
  data$y1 <- data$y1[, colSums(!is.na(data$y1)) > 0L, drop = FALSE]
  data$y0 <- data$y0[, colSums(!is.na(data$y0)) > 0L, drop = FALSE]
  
  tmp <- santosT(y1 = data$y1, y0 = data$y0, ...)
  t_ <- attr(tmp, which = 'statistic', exact = TRUE) / attr(tmp, which = 'stderr', exact = TRUE)
  
  # based on permutation
  dd <- cbind(data$y1, data$y0)
  ids <- perm_elispot(data)
  T_ <- do.call(cbind, args = lapply(ids, FUN = function(i) {
    tmp <- santosT(y1 = dd[, i, drop = FALSE], y0 = dd[, -i, drop = FALSE], ...)
    attr(tmp, which = 'statistic', exact = TRUE) / attr(tmp, which = 'stderr', exact = TRUE)
  }))
  # end of permutation
  
  ret <- maxT(t. = t_, T. = T_, ...)
  
  data$y1 <- data$y0 <- NULL
  class(data) <- 'data.frame'
  ret@design <- data
  return(ret)
}











#' @title Santos' \eqn{T}-statistic, at given time point
#' 
#' @description
#' Equation (1) and (2) of Santos' 2015 cell paper \doi{10.3390/cells4010001}.
#' 
#' @param y1,y0 \link[base]{numeric} \link[base]{matrix}-es, 
#' treatment and control responses (from an `elispot` object), respectively
#' 
#' @param null.value \link[base]{numeric} scalar \eqn{c}, 
#' as in \eqn{H_0: \bar{x}_1 - c\cdot\bar{x}_0 = 0}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @returns
#' Function [santosT] returns a \link[base]{numeric} \link[base]{vector}.
#' 
#' @importFrom DanielBiostatistics10th Gosset_Welch
#' @importFrom matrixStats rowVars
#' @importFrom stats var
#' @export
santosT <- function(
    y1, y0, 
    null.value = 1,
    ...
) {
  d1 <- dim(y1)
  m1 <- .rowMeans(y1, m = d1[1L], n = d1[2L], na.rm = TRUE)
  n1 <- .rowSums(!is.na(y1), m = d1[1L], n = d1[2L], na.rm = TRUE)
  vr1 <- rowVars(y1, na.rm = TRUE)
  # ?matrixStats::rowVars is the fastest solution I know of!!!
  
  d0 <- dim(y0)
  m0 <- .rowMeans(y0, m = d0[1L], n = d0[2L], na.rm = TRUE)
  n0 <- .rowSums(!is.na(y0), m = d0[1L], n = d0[2L], na.rm = TRUE)
  vr0 <- rowVars(y0, na.rm = TRUE)
  
  ### Equation (3) of \doi{10.3390/cells4010001} is wrong!!!! 
  ### 'the maximum standard deviation of an n = 6 binary set' # ???
  # sd_pooled <- pmax(sqrt(30)/10, sqrt( ((n1-1L)*vr1 + (n0-1L)*vr0) / (n1+n0-2L)))
  ### reason 1: need to multiply by `null.value`
  ### reason 2: should not assume equal.variance between `y1` and `y0`
  
  ret <- Gosset_Welch(v1 = vr1, v0 = vr0, c0 = null.value, n1 = n1, n0 = n0)
  # sd_pooled <- attr(ret, which = 'stderr', exact = TRUE) # this is correct!!
  attr(ret, which = 'statistic') <- (m1 - null.value * m0)
  return(ret)
}








