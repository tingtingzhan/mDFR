

#' @title Modified Distribution Free Resampling, at Two Timepoints
#' 
#' @description
#' ..
#' 
#' @param data1,data0 two `elispot` objects, 
#' at two time points \eqn{t_1} and \eqn{t_0}, respectively
#' 
#' @param ... additional parameters, such as `null.value` in function [santosT] and
#' `two.sided` for \linkS4class{maxT}
#' 
#' @examples
#' # to recreate Figure 2B
#' ds = split(santos1, f = ~ Hr + antigen)
#' maxT_santos_test(data1 = ds$`18.CEF`, data0 = ds$`0.CEF`)
#' maxT_santos_test(data1 = ds$`18.CEF`, data0 = ds$`0.CEF`, null.value = 10)
#' maxT_santos_test(data1 = ds$`18.CEF`, data0 = ds$`0.CEF`, two.sided = FALSE)
#' 
#' @importFrom methods new
#' @importFrom matrixStats colAnys
#' @export
maxT_santos_test <- function(
    data1, data0, 
    ...
) {
  
  if (nrow(data1) != nrow(data0)) {
    # timepoint1 and timepoint2 may not have same subjects!!
    d_1 <- data1; d_1$x1 <- d_1$x0 <- NULL
    d_0 <- data0; d_0$x1 <- d_0$x0 <- NULL
    tmp <- mapply(FUN = intersect, d_1, d_0)
    d_share <- as.data.frame.list(tmp[lengths(tmp) > 0L])
    data1 <- merge.data.frame(d_share, data1, by = names(d_share), all.x = TRUE)
    data0 <- merge.data.frame(d_share, data0, by = names(d_share), all.x = TRUE)
  }
  
  data1$x1 <- data1$x1[, colSums(!is.na(data1$x1)) > 0L, drop = FALSE]
  data1$x0 <- data1$x0[, colSums(!is.na(data1$x0)) > 0L, drop = FALSE]
  data0$x1 <- data0$x1[, colSums(!is.na(data0$x1)) > 0L, drop = FALSE]
  data0$x0 <- data0$x0[, colSums(!is.na(data0$x0)) > 0L, drop = FALSE]
  
  t_ <- santosT2(u = santosT(x1 = data1$x1, x0 = data1$x0, ...), 
                 v = santosT(x1 = data0$x1, x0 = data0$x0, ...))
  
  # based on permutation (remove all NA columns)
  ids1 <- perm_elispot(data1)
  ids0 <- perm_elispot(data0)
  
  fn <- function(data, id) {
    santosT(x1 = data[, id, drop = FALSE], x0 = data[, -id, drop = FALSE], ...)
  }
  tm1 <- lapply(ids1, FUN = fn, data = cbind(data1$x1, data1$x0))
  tm0 <- lapply(ids0, FUN = fn, data = cbind(data0$x1, data0$x0))
  
  ### actually [santosT2]
  t1. <- do.call(cbind, args = tm1)
  df1. <- do.call(cbind, args = lapply(tm1, FUN = attr, which = 'df', exact = TRUE))
  stderr1. <- do.call(cbind, args = lapply(tm1, FUN = attr, which = 'stderr', exact = TRUE))
  t0. <- do.call(cbind, args = tm0)
  df0. <- do.call(cbind, args = lapply(tm0, FUN = attr, which = 'df', exact = TRUE))
  stderr0. <- do.call(cbind, args = lapply(tm0, FUN = attr, which = 'stderr', exact = TRUE))
  
  id_ <- expand.grid(tm0 = seq_along(tm0), tm1 = seq_along(tm1))
  var_pooled <- ((df1.*stderr1.^2)[,id_$tm1] + (df0.*stderr0.^2)[,id_$tm0]) / (df1.[,id_$tm1] + df0.[,id_$tm0])
  T. <- (t1.[,id_$tm1] - t0.[,id_$tm0]) / sqrt(var_pooled)
  if (any(id <- colAnys(is.na(T.)))) {
    message(sprintf(fmt = '%.1f%% permutations with NA test-statistics are omitted', 1e2*mean.default(id)))
    T. <- T.[, !id]
  }
  ### but using [santosT2] is too slow on \emph{permutation-of-permutation}
  ### have to manually vectorize [santosT2] !!

  # combine `data1` and `data0` for output
  d1 <- data1; d1$x0 <- d1$x1 <- NULL
  d0 <- data0; d0$x0 <- d0$x1 <- NULL
  if (!identical(names(d1), names(d0))) stop('`elispot` at two time points must have same design')
  d <- mapply(FUN = function(c1, c0) {
    if (anyNA(c1) || anyNA(c0)) stop('does not allow NA in experiment design')
    if (all(c1 == c0)) return(c1)
    return(paste(c1, c0, sep = ' vs. '))
  }, c1 = d1, c0 = d0, SIMPLIFY = FALSE)
  tmp <- lapply(d, FUN = unique.default)
  
  ag0 <- list(...)[c('two.sided')]
  return(do.call(new, args = c(list(
    Class = 'maxT', 
    t. = t_, T. = T.,
    design = as.data.frame.list(d),
    name = paste(unlist(tmp[lengths(tmp) == 1L], use.names = FALSE), collapse = '; ')
  ), ag0[lengths(ag0, use.names = FALSE) > 0L])))
  
}










