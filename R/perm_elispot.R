
#' @title Permutation Indices of Treatment in `elispot`
#' 
#' @description
#' Permuted indices of treatment in a `elispot` object
#' 
#' @note
#' We follow the terminology used in Westfall & Young's \linkS4class{maxT} algorithm
#' (Box 2, page 82 of \doi{10.1214/ss/1056397487}), which is *permutation*.
#' 
#' In **R** convention, this is indeed a *combination* \link[utils]{combn}.
#' 
#' @param data a `elispot` object
#' 
#' @returns 
#' Function [perm_elispot] returns a \link[base]{list} of \link[base]{integer} \link[base]{vector}s.
#' 
#' @importFrom utils combn
#' @keywords internal
#' @export
perm_elispot <- function(data) {
  n1 <- dim(data$x1)[2L]
  n0 <- dim(data$x0)[2L]
  combn(n1 + n0, m = n1, simplify = FALSE)
}
