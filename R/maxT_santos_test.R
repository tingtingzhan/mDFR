

#' @title Modified Distribution Free Resampling, at Two Timepoints
#' 
#' @description
#' ..
#' 
#' @param data1,data0 two `elispot` objects, 
#' at two time points \eqn{t_1} and \eqn{t_0}, respectively
#' 
#' @param ... additional parameters, such as `null.value` in function [santosTm] and
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
  
  t_ <- santosTm(data1 = data1, data0 = data0, ...)
  
  # based on permutation (remove all NA columns)
  dd1 <- cbind(data1$y1, data1$y0)
  dd0 <- cbind(data0$y1, data0$y0)
  n1 <- length(ids1 <- perm_elispot(data1))
  n0 <- length(ids0 <- perm_elispot(data0))
  
  T_ <- list()
  for (i1 in seq_len(n1)) {
    id1 <- ids1[[i1]]
    for (i0 in seq_len(n0)) {
      id0 <- ids0[[i0]]
      new_ <- tryCatch(expr = santosTm(
        y11 = dd1[, id1, drop = FALSE], 
        y10 = dd1[, -id1, drop = FALSE], 
        y01 = dd0[, id0, drop = FALSE], 
        y00 = dd0[, -id0, drop = FALSE], 
        ...
      ), error = function(e) NULL)
      T_ <- c(T_, list(new_)) # ?base::cbind might be slow after `T_` gets big
      message('\r', i1, '/', n1, ' \u00d7 ', i0, '/', n0, ' permutation done!   ', sep = '', appendLF = FALSE)
    }
    # gc() # not needed!!
  } # slow, but not too slow at all
  message()
  # end of permutation
  
  id <- vapply(T_, FUN = anyNA, FUN.VALUE = NA)
  if (any(id)) {
    # at least one hypothesis has permuted treatment or permuted control being all-equal
    message(sprintf(fmt = '%.1f%% permutations with NA t-statistics are omitted', 1e2*mean.default(id)))
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





#' @title Santos' \eqn{T}-statistic, between time point
#' 
#' @description
#' Equation (6) and (7) of Santos' 2015 cell paper \doi{10.3390/cells4010001}.
#' 
#' @param data1,data0 (optional) two `elispot` objects, 
#' at two time points \eqn{t_1} and \eqn{t_0}, respectively
#' 
#' @param null.value \link[base]{numeric} scalar \eqn{c},
#' as in \eqn{H_0: (\bar{x}_{1,t_1} - c\cdot\bar{x}_{0,t_1}) - (\bar{x}_{1,t_0} - c\cdot\bar{x}_{0,t_0}) = 0}
#' 
#' @param y11 \link[base]{numeric} \link[base]{matrix}, time \eqn{t_1} treatment response 
#' 
#' @param y01 \link[base]{numeric} \link[base]{matrix}, time \eqn{t_0} treatment response 
#' 
#' @param y10 \link[base]{numeric} \link[base]{matrix}, time \eqn{t_1} control response 
#' 
#' @param y00 \link[base]{numeric} \link[base]{matrix}, time \eqn{t_0} control response 
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @returns
#' Function [santosTm] returns a \link[base]{numeric} \link[base]{vector}.
#' 
#' @importFrom DanielBiostatistics10th Gosset_Welch
#' @importFrom matrixStats rowVars
#' @importFrom stats var
#' @export
santosTm <- function(
    data1, data0, 
    null.value = 1,
    y11 = data1$y1, # time t1, treatment
    y01 = data0$y1, # time t0, treatment
    y10 = data1$y0, # time t1, control
    y00 = data0$y0, # time t0, control
    ...
) {
  
  # this is super time-sensitive!!!
  
  d11 <- dim(y11)
  m11 <- .rowMeans(y11, m = d11[1L], n = d11[2L], na.rm = TRUE)
  n11 <- .rowSums(!is.na(y11), m = d11[1L], n = d11[2L], na.rm = TRUE)
  vr11 <- rowVars(y11, na.rm = TRUE)
  # ?matrixStats::rowVars is the fastest solution I know of!!!
  
  d01 <- dim(y01)
  m01 <- .rowMeans(y01, m = d01[1L], n = d01[2L], na.rm = TRUE)
  n01 <- .rowSums(!is.na(y01), m = d01[1L], n = d01[2L], na.rm = TRUE)
  vr01 <- rowVars(y01, na.rm = TRUE)
  
  d10 <- dim(y10)
  m10 <- .rowMeans(y10, m = d10[1L], n = d10[2L], na.rm = TRUE)
  n10 <- .rowSums(!is.na(y10), m = d10[1L], n = d10[2L], na.rm = TRUE)
  vr10 <- rowVars(y10, na.rm = TRUE)
  
  d00 <- dim(y00)
  m00 <- .rowMeans(y00, m = d00[1L], n = d00[2L], na.rm = TRUE)
  n00 <- .rowSums(!is.na(y00), m = d00[1L], n = d00[2L], na.rm = TRUE)
  vr00 <- rowVars(y00, na.rm = TRUE)
  
  # time 1, treatment vs. control
  df1 <- Gosset_Welch(v1 = vr11, v0 = vr10, c0 = null.value, n1 = n11, n0 = n10)
  v1 <- attr(df1, which = 'stderr2', exact = TRUE)
  
  # time 0, treatment vs. control
  df0 <- Gosset_Welch(v1 = vr01, v0 = vr00, c0 = null.value, n1 = n01, n0 = n00)
  v0 <- attr(df0, which = 'stderr2', exact = TRUE)
  
  sd_pooled <- pmax(sqrt(30)/10, # copy what authors did for equation (3)
                    sqrt( (df1*v1 + df0*v0) / (df1+df0)))
  
  ((m11 - null.value * m10) - (m01 - null.value * m00)) / sd_pooled
  
}


