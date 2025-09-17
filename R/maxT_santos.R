
#' @title Modified Distribution Free Resampling
#' 
#' @description ..
#' 
#' @param data an \linkS4class{ELISpot}
#' 
#' @param ... additional parameters, such as `null.value` in function [santosT],
#' and `two.sided` for \linkS4class{maxT}
#' 
#' @references 
#' 
#' Radleigh Santos's 2015 paper \doi{10.3390/cells4010001}
#' 
#' @keywords internal
#' @export
maxT_santos <- function(data, ...) { 
  
  x1 <- data@x1
  x0 <- data@x0
  
  # before 2025-09
  #tmp <- santosT(x1 = x1, x0 = x0, ...)
  #t_ <- tmp / attr(tmp, which = 'stderr', exact = TRUE)
  t_ <- santosT(x1 = x1, x0 = x0, ...)
  
  # based on permutation
  dd <- cbind(x1, x0)
  ids <- combn_ELISpot(data)
  tmp <- lapply(ids, FUN = \(i) {
    #tmp <- santosT(x1 = dd[, i, drop = FALSE], x0 = dd[, -i, drop = FALSE], ...)
    #tmp / attr(tmp, which = 'stderr', exact = TRUE)
    santosT(x1 = dd[, i, drop = FALSE], x0 = dd[, -i, drop = FALSE], ...)
  })
  T_ <- unlist(tmp, use.names = FALSE)
  dim(T_) <- c(nrow(data@design), length(ids))
  # end of permutation
  
  tmp <- data@design |>
    lapply(FUN = \(i) {
      if (is.factor(i)) i <- as.character.factor(i)
      unique(i)
    })
  
  ag0 <- list(...)[c('two.sided')]
  return(do.call(new, args = c(list(
    Class = 'maxT', 
    t. = t_, T. = T_,
    design = data@design,
    name = paste(unlist(tmp[lengths(tmp) == 1L], use.names = FALSE), collapse = '; ')
  ), ag0[lengths(ag0, use.names = FALSE) > 0L])))
  
}


















