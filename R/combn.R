
#' @title Permuted Indices for Treatment
#' 
#' @description
#' Permuted indices for treatment.
#' 
#' @param x see **Usage**
#' 
#' @returns 
#' The `S3` generic function [combn_()] returns 
#' a \link[base]{list} of 
#' \link[base]{integer} \link[base]{vector}s,
#' as returned by the function \link[utils]{combn}.
#' 
#' @note
#' Westfall & Young's \linkS4class{maxT} algorithm
#' (Box 2, page 82 of \doi{10.1214/ss/1056397487}) 
#' used the term *permutation*.
#' 
#' However, according to **R** nomenclature, 
#' this is a *combination* 
#' (via the function \link[utils]{combn}).
#' 
#' The author uses the function name [combn_()] because,
#' unfortunately, the function \link[utils]{combn} is 
#' not an `S3` generic function.
#' 
#' @keywords internal
#' @importFrom utils combn
#' @export
combn_ <- function(x) UseMethod(generic = 'combn_')


#' @export
combn_.ELISpot <- function(x) {
  n1 <- ncol(x@x1)
  n0 <- ncol(x@x0)
  combn(x = n1 + n0, m = n1, simplify = FALSE)
}
