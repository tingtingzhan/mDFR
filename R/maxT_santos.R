
#' @title Modified Distribution Free Resampling
#' 
#' @description ..
#' 
#' @param x an \linkS4class{ELISpot}
#' 
#' @param ... additional parameters, such as `null.value` in function [santosT],
#' and `two.sided` for \linkS4class{maxT}
#' 
#' @references 
#' 
#' Radleigh Santos's 2015 paper \doi{10.3390/cells4010001}
#' 
#' @keywords internal
#' @name maxT
#' @export
maxT <- function(x, ...) UseMethod(generic = 'maxT')


#' @rdname maxT
#' @export maxT.santosT
#' @export
maxT.santosT <- function(x, ...) {
  
  x1 <- x@x1
  x0 <- x@x0
  
  dd <- cbind(x1, x0)
  ids <- combn_ELISpot(x)
  tmp <- lapply(ids, FUN = \(i) {
    santosT.matrix(x = dd[, i, drop = FALSE], x0 = dd[, -i, drop = FALSE], ...)
  })
  T_ <- unlist(tmp, use.names = FALSE)
  dim(T_) <- c(nrow(x@data@design), length(ids))
  # end of permutation
  
  tmp <- x@data@design |>
    lapply(FUN = \(i) {
      if (is.factor(i)) i <- as.character.factor(i)
      unique(i)
    })
  
  ag0 <- list(...)[c('two.sided')]
  
  return(do.call(new, args = c(list(
    Class = 'maxT', 
    t. = x@.Data, T. = T_,
    design = x@data@design,
    name = paste(unlist(tmp[lengths(tmp) == 1L], use.names = FALSE), collapse = '; ')
  ), ag0[lengths(ag0, use.names = FALSE) > 0L])))
  
}



















