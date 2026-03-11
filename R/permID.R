


#' @title Permuted Indices of Treatment
#' 
#' @description
#' Permuted indices of treatment.
#' 
#' @note
#' Westfall & Young's \linkS4class{maxT} algorithm
#' (Box 2, page 82 of \doi{10.1214/ss/1056397487}) used the term *permutation*.
#' 
#' However, in **R** convention, this is a *combination* (via the function \link[utils]{combn}).
#' 
#' @param x see **Usage**
#' 
#' @returns 
#' The `S3` generic function [permID()] returns a \link[base]{list} of \link[base]{integer} \link[base]{vector}s.
#' 
#' @keywords internal
#' @importFrom utils combn
#' @export
permID <- function(x) UseMethod(generic = 'permID')


#' @export
permID.ELISpot <- function(x) {
  n1 <- dim(x@x1)[2L]
  n0 <- dim(x@x0)[2L]
  combn(n1 + n0, m = n1, simplify = FALSE)
}
