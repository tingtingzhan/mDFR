
#' @title Permutation Indices of Treatment in `elispot`
#' 
#' @description
#' Permuted indices of treatment in a `elispot` object
#' 
#' @note
#' Westfall & Young's \linkS4class{maxT} algorithm
#' (Box 2, page 82 of \doi{10.1214/ss/1056397487}) says *permutation*.
#' 
#' In **R** convention, this is a *combination* \link[utils]{combn}.
#' 
#' @param data a `elispot` object
#' 
#' @returns 
#' Function [combn_elispot()] returns a \link[base]{list} of \link[base]{integer} \link[base]{vector}s.
#' 
#' @keywords internal
#' @importFrom utils combn
#' @export
combn_elispot <- function(data) {
  n1 <- dim(data@x1)[2L]
  n0 <- dim(data@x0)[2L]
  combn(n1 + n0, m = n1, simplify = FALSE)
}
