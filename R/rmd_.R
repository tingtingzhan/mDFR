

#' @title rmd_.maxT
#' 
#' @description
#' ..
#' 
#' @param x \linkS4class{maxT}
#' 
#' @param xnm ..
#' 
#' @param autofold \link[base]{logical} scalar
#' 
#' @param ... ..
#' 
#' @export rmd_.maxT
#' @export
rmd_.maxT <- function(x, xnm, autofold = TRUE, ...) {
  return(c(
    if (autofold) sprintf(fmt = '<details><summary>**Expand for `%s`**</summary>', x@name),
    '```{r results = \'asis\'}', 
    sprintf(fmt = 'z <- as.data.frame(%s)', xnm), # my ?mDFR::as.data.frame.maxT
    'z$adjp <- format_pval(z$adjp)',
    'z$tstat <- round(z$tstat, digits = 3L)',
    'if (length(z[[\'abs(tstat)\']])) z[[\'abs(tstat)\']] <- round(z[[\'abs(tstat)\']], digits = 3L)',
    'reactable(z)',
    '```', 
    sprintf(fmt = '```{r results = \'asis\', fig.height = %.1f, fig.width = %.1f}', 7, 6), 
    sprintf(fmt = 'suppressWarnings(grid::grid.draw(%s@gtable))', xnm), 
    # ?grid:::grid.draw.gTree; 
    # workhorse ?grid:::drawGTree
    # this step is the slowest for big data
    '```',
    '</details>',
    '<any-text>'
  ))
}


