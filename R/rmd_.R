

#' @title rmd_.maxT
#' 
#' @description
#' ..
#' 
#' @param x \linkS4class{maxT}
#' 
#' @param xnm ..
#' 
#' @param type ..
#' 
#' @param ... ..
#' 
#' @export rmd_.maxT
#' @export
rmd_.maxT <- function(x, xnm, type, ...) {
  return(c(
    if (type == 'html') sprintf(fmt = '<details><summary>**Expand for `%s`**</summary>', x@name),
    '```{r results = \'asis\'}', 
    sprintf(fmt = '.x1 <- as.data.frame(%s)', xnm), # my ?mDFR::as.data.frame.maxT
    '.x1$adjp <- format_pval(.x1$adjp)',
    '.x1$tstat <- round(.x1$tstat, digits = 3L)',
    'if (length(.x1[[\'abs(tstat)\']])) .x1[[\'abs(tstat)\']] <- round(.x1[[\'abs(tstat)\']], digits = 3L)',
    'reactable(.x1)',
    '```', 
    sprintf(fmt = '```{r results = \'asis\', fig.height = %.1f, fig.width = %.1f}', 7, 6), 
    sprintf(fmt = 'suppressWarnings(grid::grid.draw(%s@gtable))', xnm),
    '```',
    '</details>',
    '<any-text>'
  ))
}


