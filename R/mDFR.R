
#' @title Modified Distribution Free Resampling
#' 
#' @description ..
#' 
#' @param data an `elispot`
#' 
#' @param ... potential parameters
#' 
#' @examples 
#' mDFR(santos1, null.value = 1)
#' mDFR(santos1, null.value = 2)
#' 
#' @references 
#' 
#' Radleigh Santos's 2015 paper \doi{10.3390/cells4010001}
#' 
#' @export
mDFR <- function(data, ...) {
  ds <- split.elispot(data, ...)
  ret0 <- lapply(ds, FUN = function(i, ...) {
    maxT_santos(i, ...)
  }, ...)
  ret <- do.call(rbind.data.frame, args = c(ret0, list(make.row.names = FALSE)))
  # stopifnot(identical(class(ret), c('DFR', 'data.frame')))
  return(ret)
}

#' @title Modified Distribution Free Resampling, at Two Timepoints
#' 
#' @description
#' ..
#' 
#' @param data1,data0 `elispot`s at two time points \eqn{t_1} and \eqn{t_0}, respectively
#' 
#' @param null.value ..
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @export
mDFR.test <- function(
    data1, data0, 
    null.value = 1,
    ...
) {
  
  t_ <- santosTm(data1 = data1, data0 = data0, null.value = null.value)
  
  # based on permutation
  dd1 <- cbind(data1$y1, data1$y0)
  dd0 <- cbind(data0$y1, data0$y0)
  n1 <- length(ids1 <- combn_elispot(data1))
  n0 <- length(ids0 <- combn_elispot(data0))
  resamp <- NULL
  for (i1 in seq_len(n1)) {
    id1 <- ids1[[i1]]
    for (i0 in seq_len(n0)) {
      id0 <- ids0[[i0]]
      new_ <- tryCatch(expr = santosTm(
        y11 = dd1[, id1, drop = FALSE], 
        y10 = dd1[, -id1, drop = FALSE], 
        y01 = dd0[, id0, drop = FALSE], 
        y00 = dd0[, -id0, drop = FALSE], 
        null.value = null.value), error = function(e) NULL)
      resamp <- cbind(resamp, new_)
      message('\r', i1, '/', n1, ' \u00d7 ', i0, '/', n0, ' resample done!   ', sep = '', appendLF = FALSE)
    }
  } # slow, but not too slow at all
  message()
  # end of permutation
  
  # combine `data1` and `data0` for output
  d1 <- data1; d1$y0 <- d1$y1 <- NULL
  d0 <- data0; d0$y0 <- d0$y1 <- NULL
  if (!identical(names(d1), names(d0))) stop('`elispot` at two time points must have same design')
  d <- mapply(FUN = function(c1, c0) {
    if (anyNA(c1) || anyNA(c0)) stop('does not allow NA in experiment design')
    if (all(c1 == c0)) return(c1)
    return(paste(c1, c0, sep = ' vs. '))
  }, c1 = d1, c0 = d0, SIMPLIFY = FALSE)
  
  ret <- data.frame(
    d, # ?base::data.frame handles ?base::list correctly
    tstat = t_,
    adjp = maxT_(t. = t_, T. = resamp)
  )
  class(ret) <- c('DFR', class(ret))
  return(ret)
  
}




# @rdname maxT_santos
# @export 
maxT_santos <- function(
    data, 
    null.value = 1, 
    ...
) {
  
  t_ <- santosT(y1 = data$y1, y0 = data$y0, null.value = null.value)
  
  # based on permutation
    dd <- cbind(data$y1, data$y0)
    ids <- combn_elispot(data)
    resamp <- do.call(cbind, args = lapply(ids, FUN = function(i) {
      santosT(y1 = dd[, i, drop = FALSE], y0 = dd[, -i, drop = FALSE], null.value = null.value)
    }))
  # end of permutation
  
  ret <- data.frame(
    data,
    tstat = t_,
    adjp = maxT_(t. = t_, T. = resamp)
  )
  class(ret) <- c('DFR', class(ret))
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




#' @title Santos' \eqn{T}-statistic, between time point
#' 
#' @description
#' Equation (6) and (7) of Santos' 2015 cell paper \doi{10.3390/cells4010001}.
#' 
#' @param data1,data0 `elispot`s at two time points \eqn{t_1} and \eqn{t_0}, respectively
#' 
#' @param null.value \link[base]{numeric} scalar \eqn{\mu_0}
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
  df1 <- Gosset_Welch(v1 = vr11, v0 = vr10, n1 = n11, n0 = n10)
  v1 <- attr(df1, which = 'stderr2', exact = TRUE)
  
  # time 0, treatment vs. control
  df0 <- Gosset_Welch(v1 = vr01, v0 = vr00, n1 = n01, n0 = n00)
  v0 <- attr(df0, which = 'stderr2', exact = TRUE)
  
  sd_pooled <- pmax(sqrt(30)/10, # copy what authors did for equation (3)
                    sqrt( (df1*v1 + df0*v0) / (df1+df0)))
  
  ((m11 - null.value * m10) - (m01 - null.value * m00)) / sd_pooled
  
}





