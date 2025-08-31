

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
#' ds = split(santos1, f = ~ Hr + antigen)
#' list(
#'  '`maxT`' = maxT_santos_test(data1 = ds$`18.CEF`, data0 = ds$`0.CEF`)
#' ) |> rmd.tzh::render_(file = 'maxT')
#' 
#' @keywords internal
#' @importFrom methods new
#' @importFrom rmd.tzh md_
#' @importClassesFrom rmd.tzh md_lines
#' @export md_.maxT
#' @export
md_.maxT <- function(x, xnm, ...) {
  
  z1 <- c(
    # ?rmd.tzh::md_.reactable ready
    '```{r}', 
    '#| echo: false',
    sprintf(fmt = 'reactable_maxT(%s)', xnm),
    '```'
  ) |> new(Class = 'md_lines')
  
  z2 <- c(
    '```{r}',
    '#| echo: false',
    sprintf(fmt = '#| fig-height: %.1f', 7),
    sprintf(fmt = '#| fig-width: %.1f', 6),
    # sprintf(fmt = 'suppressWarnings(grid::grid.draw(%s@gtable))', xnm), # REMOVE @gtable!!!
    xnm |> sprintf(fmt = 'autoplot(%s)'),
    # ?grid:::grid.draw.gTree; 
    # workhorse ?grid:::drawGTree
    # this step is the slowest for big data
    '```'
  ) |> new(Class = 'md_lines')
  
  c(z1, z2) # ?rmd.tzh::c.md_lines
  
}



