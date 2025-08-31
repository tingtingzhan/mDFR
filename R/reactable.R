

#' @title Create \link[reactable]{reactable} from \linkS4class{maxT} Object
#' 
#' @description
#' To create a \link[reactable]{reactable} from a \linkS4class{maxT} object.
#' 
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
#' @importFrom rmd.tzh label_pvalue_sym
#' @export
reactable_maxT <- function(x, ...) {
  
  z <- as.data.frame.maxT(x)
  z$adjp <- z$adjp |> 
    label_pvalue_sym()()
  z$tstat <- round(z$tstat, digits = 3L)
  if (length(z[['abs(tstat)']])) z[['abs(tstat)']] <- round(z[['abs(tstat)']], digits = 3L)
  reactable(z, ...)

}
