

#' @title Create \link[reactable]{reactable} from \linkS4class{maxT} Object
#' 
#' @description
#' To create a \link[reactable]{reactable} from a \linkS4class{maxT} object.
#' 
#' @param x a \linkS4class{maxT} object
#' 
#' @param ... additional parameters of function \link[reactable]{reactable}
#' 
#' @returns
#' Function [reactable_maxT()] returns a \link[reactable]{reactable} object.
#' 
#' @keywords internal
#' @importFrom reactable reactable
#' @export
reactable_maxT <- function(x, ...) {
  
  x |>
    as.data.frame.maxT() |>
    reactable(...)

}
