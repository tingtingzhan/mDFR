
#' @title Modified Distribution Free Resampling
#' 
#' @description ..
#' 
#' @param data an `elispot`
#' 
#' @param ... additional parameters, such as `null.value` in function [santosT],
#' and `two.sided` for \linkS4class{maxT}
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
#' @keywords internal
#' @export
maxT_santos <- function(data, ...) { 
  
  x1 <- data@x1
  x0 <- data@x0
  
  tmp <- santosT(x1 = x1, x0 = x0, ...)
  t_ <- tmp / attr(tmp, which = 'stderr', exact = TRUE)
  
  # based on permutation
  dd <- cbind(x1, x0)
  ids <- combn_elispot(data)
  tmp <- lapply(ids, FUN = \(i) {
    tmp <- santosT(x1 = dd[, i, drop = FALSE], x0 = dd[, -i, drop = FALSE], ...)
    tmp / attr(tmp, which = 'stderr', exact = TRUE)
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


















