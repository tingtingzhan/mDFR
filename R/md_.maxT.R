

#' @title md_.maxT
#' 
#' @description
#' ..
#' 
#' @param x \linkS4class{maxT} object
#' 
#' @param xnm ..
#' 
#' @param ... ..
#' 
#' @examples
#' library(rmd.tzh)
#' ds = split(santos1, f = ~ Hr + antigen)
#' list(
#'  '`maxT`' = maxT_santos_test(data1 = ds$`18.CEF`, data0 = ds$`0.CEF`)
#' ) |> render_(file = 'maxT')
#' 
#' @keywords internal
#' @importFrom rmd.tzh md_
#' @export
md_.maxT <- function(x, xnm, ...) c(
  
  # ?rmd.tzh::md_.reactable ready
  '```{r}', 
  sprintf(fmt = 'reactable_maxT(%s)', xnm),
  '```', 
  
  '```{r}',
  sprintf(fmt = '#| fig-height: %.1f', 7),
  sprintf(fmt = '#| fig-width: %.1f', 6),
  sprintf(fmt = 'suppressWarnings(grid::grid.draw(%s@gtable))', xnm), 
  # ?grid:::grid.draw.gTree; 
  # workhorse ?grid:::drawGTree
  # this step is the slowest for big data
  '```',
  '<any-text>'
)



