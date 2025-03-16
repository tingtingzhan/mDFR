

#' @title rmd_.maxT
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
#' @keywords internal
#' @export
rmd_.maxT <- function(x, xnm, ...) c(
  
  # ?rmarkdown.tzh::rmd_.reactable ready
  '```{r}', 
  sprintf(fmt = 'reactable_maxT(%s)', xnm),
  '```', 
  
  sprintf(fmt = '```{r fig.height = %.1f, fig.width = %.1f}', 7, 6), 
  sprintf(fmt = 'suppressWarnings(grid::grid.draw(%s@gtable))', xnm), 
  # ?grid:::grid.draw.gTree; 
  # workhorse ?grid:::drawGTree
  # this step is the slowest for big data
  '```',
  '<any-text>'
)



