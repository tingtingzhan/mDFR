
#' @title Modified Distribution Free Resampling
#' 
#' @description ..
#' 
#' @param data an `elispot`
#' 
#' @param null.value ..
#' 
#' @param ... potential parameters
#' 
#' @examples 
#' # see santos_playground.R
#' 
#' @examples
#' # paper page 5, 'The ELISPOT spot results were analyzed for statistical significance at 
#' #  .. a given time point by comparing the number of spots for each of the tested peptides 
#' #  .. (triplicates) with the background controls (no peptide added, six replicates) 
#' #  .. at that time point.'
#' # paper page 8, layout of Figure 2B 
#' # these are reasons why we use `f = ~ Hr + antigen`
#' lapply(split(santos1, f = ~ Hr + antigen), FUN = maxT_santos, null.value = 1)
#' lapply(split(santos1, f = ~ Hr + antigen), FUN = maxT_santos, null.value = 2)
#' lapply(split(santos2, f = ~ Hr + antigen), FUN = maxT_santos, null.value = 1)
#' lapply(split(santos2, f = ~ Hr + antigen), FUN = maxT_santos, null.value = 2)
#' 
#' @references 
#' 
#' Radleigh Santos's 2015 paper \doi{10.3390/cells4010001}
#' 
#' @export
maxT_santos <- function(data, null.value = 1, ...) {
  data$y1 <- data$y1[, colSums(!is.na(data$y1)) > 0L, drop = FALSE]
  data$y0 <- data$y0[, colSums(!is.na(data$y0)) > 0L, drop = FALSE]
  
  t_ <- santosT(y1 = data$y1, y0 = data$y0, null.value = null.value)
  
  # based on permutation
  dd <- cbind(data$y1, data$y0)
  ids <- combn_elispot(data)
  resamp <- do.call(cbind, args = lapply(ids, FUN = function(i) {
    santosT(y1 = dd[, i, drop = FALSE], y0 = dd[, -i, drop = FALSE], null.value = null.value)
  }))
  # end of permutation
  
  ret <- maxT(t. = t_, T. = resamp, ...)
  
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
#' @param y1,y0 \link[base]{numeric} \link[base]{matrix}-es, treatment and control responses, respectively
#' 
#' @param null.value \link[base]{numeric} scalar \eqn{\mu_0}, as in \eqn{H_0: \bar{x}_1 - \mu_0\bar{x}_0 = 0}
#' 
#' @returns
#' Function [santosT] returns a \link[base]{numeric} \link[base]{vector}.
#' 
#' @importFrom stats var
#' @export
santosT <- function(y1, y0, null.value) {
  m1 <- rowMeans(y1, na.rm = TRUE)
  m0 <- rowMeans(y0, na.rm = TRUE)
  n1 <- rowSums(!is.na(y1))
  n0 <- rowSums(!is.na(y0))
  vr1 <- apply(y1, MARGIN = 1L, FUN = var, na.rm = TRUE)
  vr0 <- apply(y0, MARGIN = 1L, FUN = var, na.rm = TRUE)
  sd_pooled <- pmax(sqrt(30)/10, # 'the maximum standard deviation of an n = 6 binary set' # ???
                    sqrt( ((n1-1L)*vr1 + (n0-1L)*vr0) / (n1+n0-2L)))
  (m1 - null.value * m0) / sd_pooled
}








