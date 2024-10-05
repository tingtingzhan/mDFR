

#' @title Modified Distribution Free Resampling, at Two Timepoints
#' 
#' @description
#' ..
#' 
#' @param data1,data0 two `elispot` objects, 
#' at two time points \eqn{t_1} and \eqn{t_0}, respectively
#' 
#' @param ... additional parameters, such as `null.value` in function [santosT] and
#' `two.sided` in function [maxT]
#' 
#' @examples
#' # to recreate Figure 2B
#' ds = split(santos1, f = ~ Hr + antigen)
#' maxT_santos_test(data1 = ds$`18.CEF`, data0 = ds$`0.CEF`)
#' maxT_santos_test(data1 = ds$`18.CEF`, data0 = ds$`0.CEF`, null.value = 10)
#' maxT_santos_test(data1 = ds$`18.CEF`, data0 = ds$`0.CEF`, two.sided = FALSE)
#' 
#' @export
maxT_santos_test <- function(
    data1, data0, 
    ...
) {
  
  if (nrow(data1) != nrow(data0)) {
    # timepoint1 and timepoint2 may not have same subjects!!
    d_1 <- data1; d_1$y1 <- d_1$y0 <- NULL
    d_0 <- data0; d_0$y1 <- d_0$y0 <- NULL
    tmp <- mapply(FUN = intersect, d_1, d_0)
    d_share <- as.data.frame.list(tmp[lengths(tmp) > 0L])
    data1 <- merge.data.frame(d_share, data1, by = names(d_share), all.x = TRUE)
    data0 <- merge.data.frame(d_share, data0, by = names(d_share), all.x = TRUE)
  }
  
  data1$y1 <- data1$y1[, colSums(!is.na(data1$y1)) > 0L, drop = FALSE]
  data1$y0 <- data1$y0[, colSums(!is.na(data1$y0)) > 0L, drop = FALSE]
  data0$y1 <- data0$y1[, colSums(!is.na(data0$y1)) > 0L, drop = FALSE]
  data0$y0 <- data0$y0[, colSums(!is.na(data0$y0)) > 0L, drop = FALSE]
  
  t_ <- santosTm(tm1 = santosT(y1 = data1$y1, y0 = data1$y0, ...), 
                 tm0 = santosT(y1 = data0$y1, y0 = data0$y0, ...))
  
  # based on permutation (remove all NA columns)
  n1 <- length(ids1 <- perm_elispot(data1))
  n0 <- length(ids0 <- perm_elispot(data0))
  
  fn <- function(data, id) {
    santosT(y1 = data[, id, drop = FALSE], y0 = data[, -id, drop = FALSE], ...)
  }
  tm1 <- lapply(ids1, FUN = fn, data = cbind(data1$y1, data1$y0))
  tm0 <- lapply(ids0, FUN = fn, data = cbind(data0$y1, data0$y0))
  
  T_ <- list()
  for (i1 in seq_len(n1)) {
    for (i0 in seq_len(n0)) {
      T_ <- c(T_, list(santosTm(tm1 = tm1[[i1]], tm0 = tm0[[i0]]))) # ?base::cbind might be slow after `T_` gets big
      message('\r', i1, '/', n1, ' \u00d7 ', i0, '/', n0, ' permutation done!   ', sep = '', appendLF = FALSE)
    }
  }
  # I don't know how to message if using ?base::mapply
  # .. and even ?base::.mapply may not be much faster than this for-loop
  message()
  # end of permutation
  
  id <- vapply(T_, FUN = anyNA, FUN.VALUE = NA)
  if (any(id)) {
    # at least one hypothesis has permuted treatment or permuted control being all-equal
    message(sprintf(fmt = '%.1f%% permutations with NA test-statistics are omitted', 1e2*mean.default(id)))
    T_ <- T_[!id]
  }

  ret <- maxT(t. = t_, T. = do.call(cbind, args = T_), ...)
  
  # combine `data1` and `data0` for output
  d1 <- data1; d1$y0 <- d1$y1 <- NULL
  d0 <- data0; d0$y0 <- d0$y1 <- NULL
  if (!identical(names(d1), names(d0))) stop('`elispot` at two time points must have same design')
  d <- mapply(FUN = function(c1, c0) {
    if (anyNA(c1) || anyNA(c0)) stop('does not allow NA in experiment design')
    if (all(c1 == c0)) return(c1)
    return(paste(c1, c0, sep = ' vs. '))
  }, c1 = d1, c0 = d0, SIMPLIFY = FALSE)
  tmp <- lapply(d, FUN = unique.default)
  
  ret@design <- as.data.frame.list(d)
  ret@name <- paste(unlist(tmp[lengths(tmp) == 1L], use.names = FALSE), collapse = '; ')  
  return(ret)
  
}





#' @title Santos' Test-Statistic, between time point
#' 
#' @description
#' Equation (6) and (7) of Santos' 2015 cell paper \doi{10.3390/cells4010001}.
#' 
#' @param tm1,tm0 returns from function [santosT] at two time points \eqn{t_1} and \eqn{t_2}
#' 
#' @returns
#' Function [santosTm] returns a \link[base]{numeric} \link[base]{vector}.
#' 
#' @export
santosTm <- function(
    tm1, tm0
) {
  v1 <- attr(tm1, which = 'stderr2', exact = TRUE)
  v0 <- attr(tm0, which = 'stderr2', exact = TRUE)
  stat1 <- attr(tm1, which = 'statistic', exact = TRUE)
  stat0 <- attr(tm0, which = 'statistic', exact = TRUE)
  sd_pooled <- pmax(sqrt(30)/10, # copy what authors did for equation (3)
                    sqrt( (tm1*v1 + tm0*v0) / (tm1+tm0)))
  (stat1 - stat0) / sd_pooled
}







