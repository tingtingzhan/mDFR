

#' @title Modified Distribution Free Resampling, at Two Timepoints
#' 
#' @description
#' ..
#' 
#' @param data1,data0 two \linkS4class{ELISpot} objects, 
#' at two time points \eqn{t_1} and \eqn{t_0}, respectively
#' 
#' @param ... additional parameters, such as `null.value` in function [santosT] and
#' `two.sided` for \linkS4class{maxT}
#' 
#' @keywords internal
#' @importFrom matrixStats colAnys
#' @export
maxT_santos_test <- function(
    data1, data0, 
    ...
) {
  
  if (nrow(data1@design) != nrow(data0@design)) {
    # timepoint1 and timepoint2 may not have same subjects!!
    
    # before 2025-09
    #d_1 <- data1; d_1$x1 <- d_1$x0 <- NULL
    #d_0 <- data0; d_0$x1 <- d_0$x0 <- NULL
    #tmp <- mapply(FUN = intersect, d_1, d_0)
    #d_share <- as.data.frame.list(tmp[lengths(tmp) > 0L])
    #data1 <- merge.data.frame(d_share, data1, by = names(d_share), all.x = TRUE)
    #data0 <- merge.data.frame(d_share, data0, by = names(d_share), all.x = TRUE)
    # end of before 2025-09
    
    d_1 <- data1@design
    d_0 <- data0@design
    tmp <- mapply(FUN = intersect, d_1, d_0)
    d_share <- as.data.frame.list(tmp[lengths(tmp) > 0L])
    stop('um, tzh needs to think about how to write this beautifully..')
    data1 <- merge.data.frame(d_share, data1, by = names(d_share), all.x = TRUE)
    data0 <- merge.data.frame(d_share, data0, by = names(d_share), all.x = TRUE)
    
  }
  
  nr <- nrow(data1@design) # now `nrow(data1) == nrow(data0)`
  
  x11 <- data1@x1
  x10 <- data1@x0
  x01 <- data0@x1
  x00 <- data0@x0
  
  t_ <- santosT2(u = santosT(x1 = x11, x0 = x10, ...), 
                 v = santosT(x1 = x01, x0 = x00, ...))
  
  # based on permutation (remove all NA columns)
  ids1 <- combn_ELISpot(data1)
  n1 <- length(ids1)
  ids0 <- combn_ELISpot(data0)
  n0 <- length(ids0)
  
  fn <- function(data, id) {
    santosT(x1 = data[, id, drop = FALSE], x0 = data[, -id, drop = FALSE], ...)
  }
  tm1 <- ids1 |>
    lapply(FUN = fn, data = cbind(x11, x10))
  tm0 <- ids0 |> 
    lapply(FUN = fn, data = cbind(x01, x00))
  
  ### actually [santosT2]
  t1. <- tm1 |>
    unlist(use.names = FALSE)
  df1. <- tm1 |>
    lapply(FUN = attr, which = 'df', exact = TRUE) |>
    unlist(use.names = FALSE)
  stderr1. <- tm1 |>
    lapply(FUN = attr, which = 'stderr', exact = TRUE) |>
    unlist(use.names = FALSE)
  dim(t1.) <- dim(df1.) <- dim(stderr1.) <- c(nr, n1)
  
  t0. <- tm0 |> 
    unlist(use.names = FALSE)
  df0. <- tm0 |>
    lapply(FUN = attr, which = 'df', exact = TRUE) |>
    unlist(use.names = FALSE)
  stderr0. <- tm0 |>
    lapply(FUN = attr, which = 'stderr', exact = TRUE) |>
    unlist(use.names = FALSE)
  dim(t0.) <- dim(df0.) <- dim(stderr0.) <- c(nr, n0)
  
  id_ <- expand.grid(tm0 = seq_len(n0), tm1 = seq_len(n1))
  var_pooled <- ((df1.*stderr1.^2)[,id_$tm1] + (df0.*stderr0.^2)[,id_$tm0]) / (df1.[,id_$tm1] + df0.[,id_$tm0])
  T. <- (t1.[,id_$tm1] - t0.[,id_$tm0]) / sqrt(var_pooled)
  if (any(id <- colAnys(is.na(T.)))) {
    message(sprintf(fmt = '%.1f%% permutations with NA test-statistics are omitted', 1e2*mean.default(id)))
    T. <- T.[, !id]
  }
  ### but using [santosT2] is too slow on \emph{permutation-of-permutation}
  ### have to manually vectorize [santosT2] !!

  # combine `data1` and `data0` for output
  d1 <- data1@design
  d0 <- data0@design
  if (!identical(names(d1), names(d0))) stop('`ELISpot` at two time points must have same design')
  d <- mapply(FUN = \(c1, c0) {
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










